#!/bin/sh
export VOLK_GENERIC=1
export GR_DONT_LOAD_PREFS=1
export srcdir=/home/oper/src/phasecal_dist/gr-drf/lib
export GR_CONF_CONTROLPORT_ON=False
export PATH=/home/oper/src/phasecal_dist/gr-drf/build/lib:$PATH
export LD_LIBRARY_PATH=/home/oper/src/phasecal_dist/gr-drf/build/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$PYTHONPATH
test-drf 