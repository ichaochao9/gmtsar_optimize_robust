#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

if [[ -z "${AF_PATH:-}" ]]; then
    echo "Error: AF_PATH is not set." >&2
    echo "Set it to the ArrayFire installation directory before building." >&2
    echo >&2
    echo "Example:" >&2
    echo "  AF_PATH=/opt/arrayfire bash bash_xcorr2_cl_compiler.sh" >&2
    exit 1
fi

for command_name in gcc g++ pkg-config; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "Error: required command not found: $command_name" >&2
        exit 1
    fi
done

if ! pkg-config --exists glib-2.0; then
    echo "Error: glib-2.0 development files are not available" >&2
    exit 1
fi

if [[ ! -d "$AF_PATH/include" ]]; then
    echo "Error: ArrayFire headers not found: $AF_PATH/include" >&2
    exit 1
fi

AF_LIB_DIR=""

for candidate in "$AF_PATH/lib" "$AF_PATH/lib64"; do
    if compgen -G "$candidate/libaf.so*" >/dev/null; then
        AF_LIB_DIR="$candidate"
        break
    fi
done

if [[ -z "$AF_LIB_DIR" ]]; then
    echo "Error: ArrayFire libaf.so was not found under $AF_PATH" >&2
    exit 1
fi

gcc \
    -std=gnu99 \
    -O2 \
    -Wall \
    -Wextra \
    -I. \
    $(pkg-config --cflags glib-2.0) \
    -c xcorr2_args.c \
    -o xcorr2_args.o

gcc \
    -std=gnu99 \
    -O2 \
    -Wall \
    -Wextra \
    -I. \
    $(pkg-config --cflags glib-2.0) \
    -c prm_helper.c \
    -o prm_helper.o

g++ \
    -std=c++14 \
    -O2 \
    -Wall \
    -Wextra \
    -I. \
    -I"$AF_PATH/include" \
    $(pkg-config --cflags glib-2.0) \
    xcorr2_cl2.cpp \
    xcorr2_args.o \
    prm_helper.o \
    -L"$AF_LIB_DIR" \
    -Wl,-rpath,"$AF_LIB_DIR" \
    -laf \
    $(pkg-config --libs glib-2.0) \
    -o xcorr2_cl2

echo
echo "Build successful: $SCRIPT_DIR/xcorr2_cl2"
