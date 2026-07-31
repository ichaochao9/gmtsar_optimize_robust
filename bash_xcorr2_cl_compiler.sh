#!/usr/bin/env bash
#

#
AF_PATH="/home/wafeng/soft/insar/GMTSAR/arrayfirev3.8.3/build/install"
g++ -std=c++14 -O2 -o xcorr2_cl2 xcorr2_cl2.cpp xcorr2_args.o prm_helper.o   -H  -I$AF_PATH/include -L$AF_PATH/lib     -laf $(pkg-config --cflags --libs glib-2.0)
g++ -std=c++14 -O2 -o xcorr2_cl3 xcorr2_cl3.cpp xcorr2_args.o prm_helper.o   -H  -I$AF_PATH/include -L$AF_PATH/lib     -laf $(pkg-config --cflags --libs glib-2.0)

