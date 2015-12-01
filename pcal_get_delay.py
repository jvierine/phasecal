#!/usr/bin/env python

import h5py
import digital_rf_hdf5 as drf
import matplotlib.pyplot as plt
import numpy as n
import sys, time
import stuffr
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--dir", dest="dir", type="string", default="/data/phasecal",
                  help="Directory. (default %default)")

parser.add_option("-i", "--integrate",dest="integrate", action="store",default=100, type="int", help="Integration factor (default %default)")

parser.add_option("-b", "--baseline_time",dest="baseline_time", action="store",default=-1.0, type="float", help="Time to use as reference for cable delay, unix seconds. Default, now - 5 minutes")

parser.add_option("-0", "--t0",dest="t0", action="store",default=-1.0, type="float", help="Start time in unix seconds, default: now-5 minutes")

parser.add_option("-1", "--t1",dest="t1", action="store",default=-1.0, type="float", help="End time in unix seconds, default: now")

parser.add_option("-p", "--plot",dest="plot", action="store_true",help="plot relative time delay")

parser.add_option("-o", "--overview_plot",dest="overview_plot", action="store_true",help="plot sparse overview plot")

parser.add_option("-a", "--ascii_out",dest="ascii_out", action="store_true",help="output delays in ascii")

parser.add_option("-n", "--latest",dest="latest", action="store_true",help="Latest recorded delay")

(op, args) = parser.parse_args()

d=drf.read_hdf5(op.dir)
b0=d.get_bounds("000")
b1=d.get_bounds("001")
sample_rate = 100.0
#print(b0)
#print(b1)

t_now = time.time()
if op.baseline_time < 0.0:
    op.baseline_time = t_now - 5*60.0

if op.t0 < 0.0:
    op.t0 = t_now - 5*60.0

if op.t1 < 0.0:
    op.t1 = b0[1]/sample_rate - 2

# get baseline
try:
    z0 = n.mean(d.read_vector_c81d(long(op.baseline_time*sample_rate),op.integrate,"000"))
                            
    z1 = n.mean(d.read_vector_c81d(long(op.baseline_time*sample_rate),op.integrate,"001"))
except:
    print("Couldn't find data for determining baseline delay at %s"%(stuffr.unix2datestr(op.baseline_time)))
    exit(0)

    

# phase in 5 MHz to picoseconds ( (1/5e6)/ 1e-12)
reference_delay_ps = 200e3*n.angle(z0/z1)/2.0/n.pi
#print(""reference_delay_ps)
n_samples = long(n.floor((op.t1-op.t0)*sample_rate))

# if overview plot, calculate delay sparsely over span of data
if op.latest:
    idx0 = long(b1[1]-op.integrate)
    idx1 = long(b1[1])
    z0 = n.mean(d.read_vector_c81d(b1[1]-op.integrate,op.integrate,"000"))
    z1 = n.mean(d.read_vector_c81d(b1[1]-op.integrate,op.integrate,"001"))
    delay_ps = 200e3*n.angle(z0/z1)/2.0/n.pi - reference_delay_ps 
    print("t0 %1.3f t1 %1.3f delay %1.3f reference time %1.2f"%(idx0/sample_rate,idx1/sample_rate,delay_ps,op.baseline_time))
    exit(0)
    
if op.overview_plot:
    n_overview = 300
    z0 = n.zeros(n_overview,dtype=n.complex64)
    z1 = n.zeros(n_overview,dtype=n.complex64)
    
    t0s = n.floor(sample_rate*n.linspace(op.t0,op.t1,num=n_overview))

    tvec = t0s/sample_rate

    for ni in range(n_overview):
        try:
            z0[ni] = n.mean(d.read_vector_c81d(long(t0s[ni]),op.integrate,"000"))
            z1[ni] = n.mean(d.read_vector_c81d(long(t0s[ni]),op.integrate,"001"))
        except:
            z0[ni]=n.nan
            z1[ni]=n.nan
            print("Missing data at %s"%(stuffr.unix2datestr(t0s[ni]/sample_rate)))

else:
    # get baseline
    z0 = stuffr.decimate(d.read_vector_c81d(long(op.t0*sample_rate),n_samples,"000"),dec=op.integrate)
    
    z1 = stuffr.decimate(d.read_vector_c81d(long(op.t0*sample_rate),n_samples,"001"),dec=op.integrate)
    tvec = op.integrate*n.arange(len(z0))/sample_rate + op.t0

# phase in 5 MHz to picoseconds ( (1/5e6)/ 1e-12)
delay_ps = 200e3*n.angle(z0/z1)/2.0/n.pi - reference_delay_ps

#print(stuffr.unix2date(tvec[0]))
dates = [stuffr.unix2date(ts) for ts in tvec]

if op.plot:
    plt.plot(dates,delay_ps)
    plt.ylabel("Delay (ps)")
    plt.xlabel("Time (UTC)")
    plt.show()

if op.ascii_out:
    print("# pcal out")
    print("# reference delay at %1.2f (unix seconds), ref: %1.2f (ps)"%(op.baseline_time,reference_delay_ps))
    print("# integration %1.2f (seconds)"%(op.integrate/sample_rate))
    print("# time (unix seconds), delay (ps)")
    for ti in range(len(tvec)):
        print("%1.2f %1.2f"%(tvec[ti],delay_ps[ti]))
    


