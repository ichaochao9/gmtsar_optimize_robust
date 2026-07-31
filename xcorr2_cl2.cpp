/*
 * xcorr2_cl2.cpp -- robust ArrayFire (CUDA / OpenCL / CPU) frequency-domain
 * cross-correlation for InSAR SLC image registration (GMTSAR-GPU xcorr2).
 *
 * Based on the MIT-licensed xcorr2_cl implementation developed by
 * Hao Cui and Xianjie Zha:
 * https://github.com/cuihaoleo/gmtsar_optimize
 * Reference:
 * Cui, H., & Zha, X. (2018). Parallel Image Registration Implementations for GMTSAR Package. Seismological Research
 * Letters, 89(3), 1129-1136.
 * https://doi.org/10.1785/0220170171
 *
 * Robustness modifications are maintained in this fork:
 * https://github.com/ichaochao9/gmtsar_optimize_robust
 *
 * Usage:
 *   xcorr2_cl2 master.PRM slave.PRM output.dat [options]
 *
 * Output record format:
 *   " %d %6.3lf %d %6.3lf %6.2lf \n"
 *   (loc_x, xoff, loc_y, yoff, max_corr).
 *
 * SPDX-License-Identifier: MIT
 *
 * CHANGELOG (robustness rewrite):
 *   FIX1: Slave SLC is read with its own line width s_nx (the original code
 *         used m_nx). load_slc_rows() now uses 64-bit byte offsets, clamps
 *         start_row >= 0, derives the number of available rows from the file
 *         size, zero-fills rows past EOF, and zero-fills the buffer when the
 *         read fails. It always returns an array of exactly (n_rows x nx).
 *   FIX2: All master/slave slice windows are clamped to valid ranges so no
 *         af::seq is ever constructed with out-of-range bounds. A window that
 *         cannot be placed at all (window larger than the image) degrades to
 *         a stderr warning plus a max_corr=0.00 record instead of aborting.
 *   FIX3: No assert() in the processing loop. Sub-pixel peak indices are
 *         clamped to the valid range with a validity flag, and the whole
 *         per-patch GPU section is wrapped in try/catch (af::exception and
 *         std::exception). On failure a warning is printed and the record is
 *         still written with max_corr=0.00 so downstream line counts (e.g.
 *         fitoffset.csh) stay intact -- a bad patch never aborts the run.
 *   FIX4: NaN/Inf guards: the GMTSAR debug scaling corr *= max_corr/cmax is
 *         skipped when cmax <= 0 or non-finite; the normalized correlation
 *         score is set to 0 when denom1*denom2 is zero or non-finite; xfrac /
 *         yfrac are reset to 0 when non-finite.
 *   FIX5: Memory hygiene for thousands of patches: loop-invariant arrays
 *         (checkerboard corr_mask) are allocated once, and af::sync() +
 *         af::deviceGC() are called every 64 patches to release cached device
 *         memory (large L-band windows otherwise exhaust GPU memory).
 *   FIX6: Location sampling keeps the x_inc/y_inc scheme of xcorr2.c but
 *         every master patch center is clamped so the master window always
 *         lies inside [0, m_nx/m_ny); slave windows are clamped per FIX2.
 *   FIX7: Numerical semantics of xcorr2.c / xcorr.c are preserved:
 *         amplitude -> demean -> border mask (zero c2 borders of width
 *         xsearch/ysearch) -> 2D FFT -> checkerboard sign shift -> multiply
 *         by conj -> IFFT normalized by 1/(nx_win*ny_win) -> crop the
 *         (ysearch..ysearch+ny_corr-1, xsearch..xsearch+nx_corr-1) surface ->
 *         argmax -> normalized score
 *         100*|sum(core1*core2)|/(norm(core1)*norm(core2)) -> optional range
 *         DFT interpolation by ri -> optional sub-pixel DFT interpolation
 *         (corr2^0.25 sharpening, interp_factor) -> same output equations
 *         xoff = x_offset - (xpeak + xfrac)/ri,
 *         yoff = y_offset - (ypeak + yfrac) + loc_y*astretcha.
 */

#include <complex.h>

extern "C" {
#include "xcorr2.h"
#include "xcorr2_args.h"
}

#undef complex
#include <arrayfire.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <complex>

// Number of patches between device garbage collections (FIX5).
#define XCORR2_GC_INTERVAL 64

/*
 * load_slc_rows: read n_rows complex<int16> rows of width nx starting at row
 * `start` from an open SLC stream. Always returns an af::array of exactly
 * (n_rows x nx) complex<float> (dim0 = rows, dim1 = columns).
 *
 * FIX1: 64-bit byte offset, start_row clamped to >= 0, the number of rows
 * actually available is derived from the file size, rows past EOF (and rows
 * before row 0) are zero-filled, and the buffer is zero-filled when the read
 * itself fails -- downstream code never sees garbage and never has to deal
 * with a short/empty array.
 */
af::array load_slc_rows(std::ifstream &fin, int start, int n_rows, int nx) {
    if (n_rows <= 0 || nx <= 0)
        throw std::invalid_argument("load_slc_rows: n_rows and nx must be positive");

    // Zero-initialized interleaved (re, im) int16 buffer: missing data stays 0.
    std::vector<int16_t> buf((size_t)n_rows * (size_t)nx * 2, 0);

    const int64_t row_bytes = (int64_t)nx * 2 * (int64_t)sizeof(int16_t);
    // FIX1: clamp start_row >= 0; leading rows of the window stay zero-filled.
    const int start_row = start < 0 ? 0 : start;
    const int skipped_rows = start_row - start;

    if (fin) {
        // Determine the file size (64-bit) to know how many rows exist.
        fin.clear();
        fin.seekg(0, std::ios::end);
        const int64_t file_size = (int64_t)fin.tellg();

        // FIX1: 64-bit byte offset -- the original used a 32-bit long
        // expression that overflows for large L-band SLCs.
        const int64_t offset = (int64_t)nx * start_row * 2 * (int64_t)sizeof(int16_t);

        if (file_size > offset) {
            const int64_t avail_rows = (file_size - offset) / row_bytes;
            const int want_rows = n_rows - skipped_rows;
            const int rows_to_read = (int)std::min<int64_t>(avail_rows, want_rows);

            if (rows_to_read > 0) {
                fin.seekg(offset, std::ios::beg);
                int16_t *dst = buf.data() + (size_t)skipped_rows * (size_t)nx * 2;
                fin.read((char *)dst, (std::streamsize)rows_to_read * row_bytes);
                if (!fin) {
                    // FIX1: check stream state after read; zero-fill on failure.
                    std::fill(dst, dst + (size_t)rows_to_read * (size_t)nx * 2, (int16_t)0);
                    fin.clear();
                }
            }
        }
    }

    af::array af_buf(2, nx, n_rows, buf.data());
    af::array dest = af::complex(af_buf(0, af::span, af::span), af_buf(1, af::span, af::span));
    dest = af::moddims(dest, nx, n_rows);

    return af::transpose(dest); // -> (n_rows x nx), complex<float>
}

/*
 * dft_interpolate: zero-padding DFT interpolation by (scale_h, scale_w).
 * Unchanged from the original implementation (FIX7 semantics); all af::seq
 * bounds below are provably in-range for height,width >= 1:
 *   left  [0, width/2], right [width/2, width-1],
 *   up    [0, height/2-1], down [height/2, height-1].
 */
af::array dft_interpolate(const af::array &in, int scale_h, int scale_w) {
    int height = in.dims(0);
    int width = in.dims(1);
    int out_height = height * scale_h;
    int out_width = width * scale_w;

    af::array in_fft = af::dft(in);
    af::array out_fft = af::constant(af::cdouble(0, 0), out_height, out_width, c32);

    af::seq left = af::seq(0, width/2);
    af::seq &out_left = left;
    af::seq right = af::seq(width/2, width-1);
    af::seq out_right = af::seq(out_width-width/2, out_width-1);
    af::seq up = af::seq(0, height/2-1);
    af::seq &out_up = up;
    af::seq down = af::seq(height/2, height-1);
    af::seq out_down = af::seq(out_height-height/2, out_height-1);

    in_fft(af::span, width/2) /= 2.0;
    out_fft(out_up, out_left) = in_fft(up, left);
    out_fft(out_up, out_right) = in_fft(up, right);
    out_fft(out_down, out_left) = in_fft(down, left);
    out_fft(out_down, out_right) = in_fft(down, right);

    return af::idft(out_fft, 1.0/(height * width), out_fft.dims());
}

void arrayfire_init(const st_xcorr_args &args) {
    switch (args.device) {
        case XCORR2_DEVICE_CUDA:
            af::setBackend(AF_BACKEND_CUDA);
            break;
        case XCORR2_DEVICE_OPENCL:
            af::setBackend(AF_BACKEND_OPENCL);
            break;
        case XCORR2_DEVICE_CPU:
            af::setBackend(AF_BACKEND_CPU);
            break;
        default:
            af::setBackend(AF_BACKEND_DEFAULT);
    }
    af::info();
}

int main(int argc, char **argv) {
    st_xcorr_args args;
    st_xcorr xcorr;

    parse_opts(&args, argc, argv);
    apply_args(&args, &xcorr);
    arrayfire_init(args);

    const int xsearch = xcorr.xsearch;
    const int ysearch = xcorr.ysearch;
    const int nx_corr = xcorr.xsearch * 2;
    const int nx_win = nx_corr * 2;
    const int ny_corr = xcorr.ysearch * 2;
    const int ny_win = ny_corr * 2;

    // Startup sanity checks (config errors, reported before processing begins).
    if (xsearch < 1 || ysearch < 1) {
        fprintf(stderr, "xcorr2_cl: xsearch/ysearch must be >= 1\n");
        return 1;
    }
    // x_inc/y_inc below divide by (nxl+3) and (nyl+1); reject configs where
    // the patch counts are invalid or the divisors would be zero.
    if (xcorr.nxl < 1 || xcorr.nyl < 1 ||
        xcorr.nxl + 3 == 0 || xcorr.nyl + 1 == 0) {
        fprintf(stderr,
                "xcorr2_cl2: invalid nxl/nyl (%d, %d): sampling divisor would be zero\n",
                xcorr.nxl, xcorr.nyl);
        return 1;
    }
    // ri is a divisor in the output equation; guard against division by zero.
    const float ri = xcorr.ri > 0 ? (float)xcorr.ri : 1.0f;

    // FIX6: keep the x_inc/y_inc sampling scheme of xcorr2.c.
    const int x_inc = (xcorr.m_nx - 2*(xcorr.xsearch + nx_corr)) / (xcorr.nxl + 3);
    const int y_inc = (xcorr.m_ny - 2*(xcorr.ysearch + ny_corr)) / (xcorr.nyl + 1);

    std::ifstream f1(xcorr.m_path, std::ios::binary);
    if (!f1) {
        fprintf(stderr, "xcorr2_cl: failed to open master SLC image %s\n", xcorr.m_path);
        return 1;
    }
    std::ifstream f2(xcorr.s_path, std::ios::binary);
    if (!f2) {
        fprintf(stderr, "xcorr2_cl: failed to open slave SLC image %s\n", xcorr.s_path);
        return 1;
    }

    FILE *fout = fopen(xcorr.o_name, "w");
    if (!fout) {
        perror("xcorr2_cl2: failed to open out file");
        return 1;
    }

    // FIX5: loop-invariant checkerboard sign mask, allocated once.
    // Built exactly as in the original code: sign = ((i+j)&1) ? -1 : 1 over
    // the ny_win x nx_win grid (i indexes columns, j indexes rows of the
    // column-major buffer).
    std::vector<int> corr_mask_arr((size_t)nx_win * (size_t)ny_win);
    for (int i=0; i<nx_win; i++)
        for (int j=0; j<ny_win; j++)
            corr_mask_arr[(size_t)i*ny_win + j] = ((i + j) & 1) ? -1 : 1;
    af::array corr_mask(ny_win, nx_win, corr_mask_arr.data());
    corr_mask_arr.clear();

    // FIX2: a window that can never be placed (window larger than the image)
    // must not abort the run; warn once and emit max_corr=0.00 records.
    const bool master_placeable = (xcorr.m_nx >= nx_win) && (xcorr.m_ny >= ny_win);
    const bool slave_placeable  = (xcorr.s_nx >= nx_win);
    bool warned_unplaceable = false;
    bool warned_rowload = false;

    long patch_count = 0;

    for (int j=1; j<=xcorr.nyl; j++) {
        int loc_y = ny_win + j * y_inc;
        // FIX6: guarantee the master window rows stay inside [0, m_ny).
        if (master_placeable) {
            if (loc_y < ny_win/2) loc_y = ny_win/2;
            else if (loc_y > xcorr.m_ny - ny_win/2) loc_y = xcorr.m_ny - ny_win/2;
        }
        const int slave_loc_y = (int)((1+xcorr.astretcha)*loc_y + xcorr.y_offset);

        // Load row blocks. FIX1: the slave block is read with s_nx, not m_nx.
        // FIX1/FIX2: load_slc_rows clamps the start row and zero-fills rows
        // beyond EOF, so slave_loc_y - ny_win/2 may fall outside the image
        // without breaking anything (zeros behave like the masked borders).
        af::array m_rows, s_rows;
        bool rows_ok = true;
        try {
            m_rows = load_slc_rows(f1, loc_y - ny_win/2, ny_win, xcorr.m_nx);
            s_rows = load_slc_rows(f2, slave_loc_y - ny_win/2, ny_win, xcorr.s_nx);
        } catch (const af::exception &e) {
            // FIX3: a failed row block degrades to zero-score records, not abort.
            fprintf(stderr, "xcorr2_cl2: warning: failed to load row block at loc_y=%d: %s\n",
                    loc_y, e.what());
            rows_ok = false;
        } catch (const std::exception &e) {
            fprintf(stderr, "xcorr2_cl2: warning: failed to load row block at loc_y=%d: %s\n",
                    loc_y, e.what());
            rows_ok = false;
        }

        for (int i=2; i<=xcorr.nxl+1; i++) {
            patch_count++;

            int loc_x = nx_win + i * x_inc;
            // FIX6: guarantee the master window columns stay inside [0, m_nx).
            if (master_placeable) {
                if (loc_x < nx_win/2) loc_x = nx_win/2;
                else if (loc_x > xcorr.m_nx - nx_win/2) loc_x = xcorr.m_nx - nx_win/2;
            }
            const int slave_loc_x = (int)((1+xcorr.astretcha)*loc_x + xcorr.x_offset);

            // Best-effort output values (used for skipped/failed patches so
            // that downstream line counts stay consistent).
            double xoff = xcorr.x_offset;
            double yoff = xcorr.y_offset + loc_y * xcorr.astretcha;
            float max_corr = 0.0f;

            // FIX2: skip (with a record, not an abort) if a window can never
            // be placed inside its image.
            if (!rows_ok || !master_placeable || !slave_placeable) {
                if (!rows_ok) {
                    // REVIEWER-2: the row-block load error was already reported
                    // when the block was read; do not emit the misleading
                    // "window cannot be placed" warning for this cause.
                    if (!warned_rowload) {
                        warned_rowload = true;
                        fprintf(stderr,
                                "xcorr2_cl: warning: skipping patch: row-block load failed "
                                "(see earlier warning); writing max_corr=0.00 records\n");
                    }
                } else if (!warned_unplaceable) {
                    warned_unplaceable = true;
                    fprintf(stderr,
                            "xcorr2_cl: warning: correlation window (%d x %d) cannot be fully "
                            "placed in master (%d x %d) or slave (%d wide) image; "
                            "writing max_corr=0.00 records for affected patches\n",
                            nx_win, ny_win, xcorr.m_nx, xcorr.m_ny, xcorr.s_nx);
                }
                fprintf(fout, " %d %6.3lf %d %6.3lf %6.2lf \n",
                        loc_x, xoff, loc_y, yoff, (double)max_corr);
                if (patch_count % XCORR2_GC_INTERVAL == 0) {
                    // REVIEWER-6: GC is best-effort; a failure here must
                    // never terminate the run.
                    try { af::sync(); af::deviceGC(); } catch (...) { /* ignore */ }
                }
                continue;
            }

            // FIX3: wrap the whole per-patch GPU section in try/catch so one
            // bad patch (low-coherence L-band patches are common) never
            // terminates the run.
            try {
                // FIX2: clamp the slave window so [sx0, sx0+nx_win-1] lies in
                // [0, s_nx-1]; slave_placeable guarantees s_nx-nx_win >= 0.
                int sx0 = slave_loc_x - nx_win/2;
                if (sx0 < 0) sx0 = 0;
                else if (sx0 > xcorr.s_nx - nx_win) sx0 = xcorr.s_nx - nx_win;

                // Master window start; in-bounds by the FIX6 clamp above.
                const int mx0 = loc_x - nx_win/2;

                // All af::seq bounds below are provably valid:
                //   slice_x       [0, m_nx-1], slave_slice_x [0, s_nx-1].
                const af::seq slice_x(mx0, mx0 + nx_win - 1);
                const af::seq slave_slice_x(sx0, sx0 + nx_win - 1);

                af::array c1 = m_rows(af::span, slice_x);
                af::array c2 = s_rows(af::span, slave_slice_x);

                // FIX7: optional range DFT interpolation by factor ri, then
                // slice the center nx_win columns of the ri*nx_win wide rows.
                if (xcorr.ri > 1) {
                    const int interp_width = xcorr.ri * nx_win;

                    af::array interp1 = dft_interpolate(c1, 1, xcorr.ri);
                    af::array interp2 = dft_interpolate(c2, 1, xcorr.ri);

                    // In-bounds: nx_win*(ri-1)/2 .. nx_win*(ri+1)/2 - 1
                    // within [0, ri*nx_win - 1] for ri >= 1.
                    const af::seq x_seq(interp_width/2 - nx_win/2,
                                        interp_width/2 + nx_win/2 - 1);
                    c1 = interp1(af::span, x_seq);
                    c2 = interp2(af::span, x_seq);
                }

                // FIX7: amplitude -> demean, as in xcorr2.c.
                af::array c1r = af::abs(c1);
                af::array c2ro = af::abs(c2);

                float m1 = af::mean<float>(c1r);
                float m2 = af::mean<float>(c2ro);

                c1r -= m1;
                c2ro -= m2;

                // Border mask: zero the c2 borders of width xsearch/ysearch.
                // roi bounds [ysearch, ny_win-ysearch-1] x [xsearch,
                // nx_win-xsearch-1] are in-range because nx_win = 4*xsearch
                // and ny_win = 4*ysearch.
                af::array c2r = af::constant(0, ny_win, nx_win, f32);
                const af::seq roi_y(ysearch, ny_win - ysearch - 1);
                const af::seq roi_x(xsearch, nx_win - xsearch - 1);
                c2r(roi_y, roi_x) = c2ro(roi_y, roi_x);

                // FIX7: 2D FFT -> checkerboard sign shift -> multiply by
                // conj -> IFFT normalized by 1/(nx_win*ny_win).
                af::array c1r_fft = af::dft(c1r);
                af::array c2r_fft = af::dft(c2r);
                af::array c3r_fft = c1r_fft * corr_mask * af::conjg(c2r_fft);
                af::array c3r = af::idft(c3r_fft, 1.0/(nx_win*ny_win), c3r_fft.dims());

                // Crop the correlation surface; bounds [ysearch,
                // ysearch+ny_corr-1] x [xsearch, xsearch+nx_corr-1] are within
                // the (ny_win x nx_win) surface since ny_win = 2*ny_corr.
                af::array corr = af::abs(c3r(
                        af::seq(ysearch, ysearch + ny_corr - 1),
                        af::seq(xsearch, xsearch + nx_corr - 1)));

                // af::max on an empty array is invalid; cannot happen for
                // xsearch,ysearch >= 1, but guard anyway.
                if (corr.elements() == 0)
                    throw std::runtime_error("empty correlation surface");

                // Column-major linear index decomposition (original
                // convention): row = idx % dim0 (y), col = idx / dim0 (x).
                unsigned max_idx = 0;
                float cmax = 0.0f;
                af::max<float>(&cmax, &max_idx, corr);

                int xpeak = (int)(max_idx / (unsigned)ny_corr) - xsearch;
                int ypeak = (int)(max_idx % (unsigned)ny_corr) - ysearch;

                // Normalized correlation score over the matching core windows
                // (FIX7). Bounds: xpeak in [-xsearch, xsearch-1] so
                // xsearch+xpeak in [0, nx_corr-1] and the seq end
                // xsearch+xpeak+nx_corr-1 <= 2*nx_corr-2 < nx_win; same for y.
                af::array core1 = c1r(
                    af::seq(ysearch + ypeak, ysearch + ypeak + ny_corr - 1),
                    af::seq(xsearch + xpeak, xsearch + xpeak + nx_corr - 1));
                af::array core2 = c2r(
                    af::seq(ysearch, ysearch + ny_corr - 1),
                    af::seq(xsearch, xsearch + nx_corr - 1));

                // ArrayFire 3.x: the norm-type enumerators live at global
                // scope (namespace af only has typedefs), and af::norm
                // returns double.
                double denom1 = af::norm(core1, AF_NORM_EUCLID); // L2/Frobenius
                double denom2 = af::norm(core2, AF_NORM_EUCLID);
                double num = af::sum<float>(core1 * core2);

                // FIX4: guard zero / non-finite denominators.
                const double denom = denom1 * denom2;
                if (denom > 0.0 && std::isfinite(denom) && std::isfinite((double)num))
                    max_corr = (float)(100.0 * fabs((double)num / denom));
                else
                    max_corr = 0.0f;

                // Optional sub-pixel DFT interpolation (FIX7).
                float xfrac = 0.0f, yfrac = 0.0f;
                if (xcorr.interp_factor > 1) {
                    const int factor = xcorr.interp_factor;
                    const int nx_corr2 = xcorr.n2x;
                    const int ny_corr2 = xcorr.n2y;

                    // Boundary note: when n2 == 2*search exactly (i.e.
                    // nx_corr2 == nx_corr), the peak clamp below has an empty
                    // valid range and can produce an af::seq starting at -1,
                    // so interpolation must be skipped at equality too.
                    if (nx_corr2 < 2 || ny_corr2 < 2 ||
                        nx_corr2 >= nx_corr || ny_corr2 >= ny_corr) {
                        // Invalid sub-pixel window for this surface; keep the
                        // integer peak instead of crashing (warn once per run).
                        static bool warned_interp = false;
                        if (!warned_interp) {
                            warned_interp = true;
                            fprintf(stderr,
                                    "xcorr2_cl: warning: sub-pixel window (%d x %d) incompatible "
                                    "with correlation surface (%d x %d); skipping interpolation\n",
                                    nx_corr2, ny_corr2, nx_corr, ny_corr);
                        }
                    } else {
                        // Scale to match GMTSAR (kept from the original code).
                        // FIX4: skip the scaling when cmax <= 0 or non-finite.
                        if (cmax > 0.0f && std::isfinite((double)cmax) &&
                            std::isfinite((double)max_corr))
                            corr *= max_corr / cmax;

                        // Keep the original GMTSAR memory-violation workaround:
                        // offset the peak so the n2x/n2y sub-window stays fully
                        // inside the corr surface. After this,
                        // ypeak+ysearch in [ny_corr2/2, ny_corr-ny_corr2/2-1];
                        // the guard above guarantees ny_corr2 < ny_corr, so
                        // this clamp range is non-empty (same for x) and the
                        // slice below is provably in-bounds.
                        if (ypeak + ysearch < ny_corr2/2)
                            ypeak = ny_corr2 / 2 - ysearch;
                        else if (ypeak + ysearch >= ny_corr - ny_corr2/2)
                            ypeak = ny_corr - ny_corr2/2 - ysearch - 1;

                        if (xpeak + xsearch < nx_corr2/2)
                            xpeak = nx_corr2 / 2 - xsearch;
                        else if (xpeak + xsearch >= nx_corr - nx_corr2/2)
                            xpeak = nx_corr - nx_corr2/2 - xsearch - 1;

                        af::array corr2 = corr(
                                af::seq(ypeak + ysearch - ny_corr2/2,
                                        ypeak + ysearch + ny_corr2/2 - 1),
                                af::seq(xpeak + xsearch - nx_corr2/2,
                                        xpeak + xsearch + nx_corr2/2 - 1));
                        corr2 = af::pow(corr2, 0.25); // ^0.25 sharpening (FIX7)

                        af::array hi_corr = af::abs(dft_interpolate(corr2, factor, factor));

                        const int ny_hi = ny_corr2 * factor;
                        const int nx_hi = nx_corr2 * factor;

                        if (hi_corr.elements() == 0)
                            throw std::runtime_error("empty interpolated surface");

                        unsigned max_idx2 = 0;
                        float cmax2 = 0.0f;
                        af::max<float>(&cmax2, &max_idx2, hi_corr);

                        // REVIEWER-5: a degenerate hi-res surface (all-zero,
                        // or a non-finite maximum) yields a bogus peak2; keep
                        // the integer-peak result (xfrac = yfrac = 0) instead.
                        const bool hires_valid =
                            (cmax2 > 0.0f) && std::isfinite((double)cmax2);

                        int xpeak2 = (int)(max_idx2 / (unsigned)ny_hi) - nx_hi / 2;
                        int ypeak2 = (int)(max_idx2 % (unsigned)ny_hi) - ny_hi / 2;

                        // FIX3: the original asserted the peak bounds; clamp
                        // instead and flag invalid peaks (fracs stay 0).
                        bool peak_valid = true;
                        if (xpeak2 < -nx_hi/2) { xpeak2 = -nx_hi/2; peak_valid = false; }
                        else if (xpeak2 >= nx_hi/2) { xpeak2 = nx_hi/2 - 1; peak_valid = false; }
                        if (ypeak2 < -ny_hi/2) { ypeak2 = -ny_hi/2; peak_valid = false; }
                        else if (ypeak2 >= ny_hi/2) { ypeak2 = ny_hi/2 - 1; peak_valid = false; }

                        if (peak_valid && hires_valid) {
                            xfrac = xpeak2 / (float)factor;
                            yfrac = ypeak2 / (float)factor;
                        }
                        // FIX4: reset non-finite sub-pixel fractions to 0.
                        if (!std::isfinite((double)xfrac)) xfrac = 0.0f;
                        if (!std::isfinite((double)yfrac)) yfrac = 0.0f;
                    }
                }

                // FIX7: same output equations as xcorr2.c.
                xoff = xcorr.x_offset - ((xpeak + xfrac) / ri);
                yoff = xcorr.y_offset - (ypeak + yfrac) + loc_y * xcorr.astretcha;
            } catch (const af::exception &e) {
                // FIX3: warn and keep the best-effort zero-score record.
                fprintf(stderr,
                        "xcorr2_cl2: warning: GPU failure at patch (loc_x=%d, loc_y=%d): %s; "
                        "writing max_corr=0.00 record\n", loc_x, loc_y, e.what());
                xoff = xcorr.x_offset;
                yoff = xcorr.y_offset + loc_y * xcorr.astretcha;
                max_corr = 0.0f;
            } catch (const std::exception &e) {
                fprintf(stderr,
                        "xcorr2_cl2: warning: failure at patch (loc_x=%d, loc_y=%d): %s; "
                        "writing max_corr=0.00 record\n", loc_x, loc_y, e.what());
                xoff = xcorr.x_offset;
                yoff = xcorr.y_offset + loc_y * xcorr.astretcha;
                max_corr = 0.0f;
            }

            fprintf(fout, " %d %6.3lf %d %6.3lf %6.2lf \n",
                    loc_x, xoff, loc_y, yoff, (double)max_corr);

            // FIX5: periodically release cached device memory so thousands of
            // large L-band patches do not exhaust GPU memory.
            if (patch_count % XCORR2_GC_INTERVAL == 0) {
                // REVIEWER-6: GC is best-effort; a failure here must
                // never terminate the run.
                try { af::sync(); af::deviceGC(); } catch (...) { /* ignore */ }
            }
        }
    }

    f1.close();
    f2.close();
    fclose(fout);

    return 0;
}
