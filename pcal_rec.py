#!/usr/bin/env python
#
# Dual channel double precision digital down conversion
# record continuously to disk with digital_rf
#
# (c) 2015 Juha Vierinen
#
from PyQt4 import Qt
from gnuradio import qtgui
import sys, sip, time

from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from optparse import OptionParser

from gnuradio import filter
from gnuradio import blocks
from gnuradio.filter import firdes

import sys, time, os, math, re

import drf
import sampler_util

import numpy
import scipy
import scipy.signal

import serial

def create_simple_scope(samp_rate):
    snk = qtgui.time_sink_c(1024,samp_rate,"Phase",1)
    snk_win = sip.wrapinstance(snk.pyqwidget(), Qt.QWidget)
    snk_win.show()
    return(snk)

def create_simple_scope_f(samp_rate):
    snk = qtgui.time_sink_f(1024,samp_rate,"Phase",1)
    snk.set_y_axis(-4,4)
    snk_win = sip.wrapinstance(snk.pyqwidget(), Qt.QWidget)
    snk_win.show()
    return(snk)

debug = True

parser = OptionParser(option_class=eng_option)

parser.add_option("-m", "--mainboard", dest="main_board", type="string", default="192.168.10.2",
                  help="Mainboard address. (default %default)")

parser.add_option("-0", "--centerfreq0",dest="centerfreq0", action="store", type="string",
                  default="5e6",help="Center frequency (default %default)")

parser.add_option("-1", "--centerfreq1",dest="centerfreq1", action="store", type="string",
                  default="214791314.0*100e6/float(2**32)",help="Center frequency (default %default)")

parser.add_option("-r", "--samplerate",dest="samplerate", action="store", type="long", default=25000000,
                  help="Sample rate (default %default)")

parser.add_option("-i", "--window_len",dest="window_len", action="store", type="int", default=250000,
                  help="Averaging window length (default %default)")


parser.add_option("-n", "--name",dest="name", action="store", default="data",
                  type="string",
                  help="Prefix of data directory.")


(op, args) = parser.parse_args()
op.centerfreq0 = eval(op.centerfreq0)
op.centerfreq1 = eval(op.centerfreq1)

# create usrp source block
u = uhd.usrp_source(
   device_addr="addr=%s,recv_buff_size=1000000000"%(op.main_board),
   stream_args=uhd.stream_args(
      cpu_format="sc16",
      otw_format="sc16",
      channels=range(1),
      ),
   )

# tbd fix hardcoded directory
os.system("mkdir -p /data/phasecal/000")
os.system("mkdir -p /data/phasecal/001")

u.set_clock_source("external", 0)
u.set_time_source("external", 0)

u.set_subdev_spec("A:AB", 0)
u.set_samp_rate(op.samplerate)

u.set_auto_dc_offset(False,0)


tt = time.time()
while tt-math.floor(tt) < 0.2 or tt-math.floor(tt) > 0.3:
    tt = time.time()
    time.sleep(0.01)
print("Latching at "+str(tt))
print(tt)
u.set_time_unknown_pps(uhd.time_spec(math.ceil(tt)+1.0))

op.starttime = math.ceil(time.time())+5
op.starttime = sampler_util.find_next(op.starttime,1)

u.set_start_time(uhd.time_spec(op.starttime))

fg = gr.top_block()
dst = []

win_len = op.window_len
window = numpy.array(numpy.real(scipy.signal.blackmanharris(win_len)),dtype=numpy.float64)
window.tofile("/data/phasecal/coeffs.bin")

# two channel decimating all double precision digital down converter
# with user defined filter coefficients.
# - 100 samples per file.
# - let the block know that the sample rate is 25 MHz also
# - expects real valued samples with channel 0 and 1 as interleaved little endian 16 bit ints.
# - dies if packet drop detected.
dddc = drf.dddc("/data/phasecal/coeffs.bin", win_len, op.centerfreq0, op.centerfreq1, 100, 25000000.0)
fg.connect((u, 0), (dddc,0))

if debug:
   print "Starting at %d"%(op.starttime)
fg.start()

# die if loss of lock.
lt0 = time.time()


while(True):
    tnow = time.time()

    # if not locked to reference, die!
    if tnow - lt0 > 10:
        lock = str(u.get_mboard_sensor("ref_locked"))
        if lock != "Ref: locked":
            print(lock)
            print("Reference unlocked!")
            exit(0)
        lt0=tnow
    
    time.sleep(0.1)
