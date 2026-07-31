# Robust ArrayFire xcorr2 for GMTSAR (`xcorr2_cl2`)

This repository is an unofficial research fork of
[cuihaoleo/gmtsar_optimize](https://github.com/cuihaoleo/gmtsar_optimize).

It provides `xcorr2_cl2`, a robustness-oriented ArrayFire implementation of
frequency-domain cross-correlation for GMTSAR SLC image registration.

> Use `xcorr2_cl2` for the robust implementation maintained by this fork.
> The original upstream source files are retained for attribution, history,
> and supporting code.

## Main improvements

Compared with the original ArrayFire implementation, `xcorr2_cl2` adds:

- safer master and slave SLC input handling;
- separate slave-image dimensions and 64-bit file offsets;
- bounds checking for master and slave image windows;
- zero filling for unavailable or unreadable SLC rows;
- numerical guards for zero, NaN, and infinite values;
- per-patch exception recovery instead of aborting the entire run;
- periodic ArrayFire synchronization and device-memory cleanup.

Detailed implementation notes are included at the beginning of
[`xcorr2_cl2.cpp`](xcorr2_cl2.cpp).

## Status

This code is currently provided as a research preview.

The build script targets Linux systems with an ArrayFire installation
containing `libaf.so`. It has not yet been systematically validated with every
SAR sensor, ArrayFire version, or CUDA/OpenCL device.

Users should compare the resulting offsets with a trusted GMTSAR workflow
before using the program for production processing.

No precompiled binary or test SLC data are included.

## Requirements

- Linux
- GCC and G++
- `pkg-config`
- GLib 2 development files
- ArrayFire with at least one available backend:
  - CUDA
  - OpenCL
  - CPU
- GMTSAR-compatible PRM and complex SLC files

The original development environment used ArrayFire 3.8.3. Compatibility with
other ArrayFire 3.x versions has not yet been systematically tested.

On Debian or Ubuntu, the basic compiler and GLib dependencies can be installed
with:

```bash
sudo apt update
sudo apt install build-essential pkg-config libglib2.0-dev
```

ArrayFire and the required CUDA or OpenCL runtime must be installed separately.

## Build

Download this repository or clone it:

```bash
git clone https://github.com/ichaochao9/gmtsar_optimize_robust.git
cd gmtsar_optimize_robust
```

Set `AF_PATH` to the ArrayFire installation prefix. The directory must contain:

```text
$AF_PATH/include
$AF_PATH/lib/libaf.so
```

Alternatively, the shared library may be located under `$AF_PATH/lib64`.

Build `xcorr2_cl2` with:

```bash
AF_PATH=/opt/arrayfire bash bash_xcorr2_cl_compiler.sh
```

Replace `/opt/arrayfire` with the actual ArrayFire installation path on your
computer.

A successful build creates:

```text
xcorr2_cl2
```

in the repository directory.

## Usage

The output filename is a required third positional argument:

```bash
./xcorr2_cl2 master.PRM slave.PRM output.dat [options]
```

Example using the CUDA backend:

```bash
./xcorr2_cl2 master.PRM slave.PRM output.dat \
  -af cuda \
  -xsearch 64 \
  -ysearch 64 \
  -nx 32 \
  -ny 64 \
  -range_interp 2 \
  -interp 16
```

Other ArrayFire backends may be selected with:

```bash
-af opencl
```

or:

```bash
-af cpu
```

Display the program help with:

```bash
./xcorr2_cl2 -help
```

## Output

Each successful or recovered patch writes one ASCII record containing:

```text
loc_x  xoff  loc_y  yoff  max_corr
```

The source-code output format is:

```text
%d %6.3lf %d %6.3lf %6.2lf
```

The output can subsequently be processed using the appropriate GMTSAR offset
fitting workflow, such as `fitoffset.csh`.

## Known limitation

The region-selection options `-x0`, `-x1`, `-y0`, and `-y1` are not implemented
in this preview release. They are explicitly rejected instead of being silently
ignored.

## Attribution and citation

This work builds on the original MIT-licensed parallel registration
implementation developed by Hao Cui and Xianjie Zha:

> Cui, H., & Zha, X. (2018). Parallel Image Registration Implementations for
> GMTSAR Package. *Seismological Research Letters*, 89(3), 1129–1136.
> https://doi.org/10.1785/0220170171

Original project:

- https://github.com/cuihaoleo/gmtsar_optimize

This repository is an independent research fork and is not an official GMTSAR
release.

## License

This repository retains the original MIT License and copyright notice.
See [`LICENSE`](LICENSE) for details.
