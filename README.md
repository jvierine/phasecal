# phasecal 

This software package allows a relative cable delay to be measured with an accuracy of ~ 1 picosecond by tracking the phase of an upgoing 5 MHz signal and a downgoing 5 MHz + 1 kHz signal. The software requires a custom FPGA image for the USRP N200, which simply streams dual channel decimated raw 100 MHz ADC samples. 

Hardware
--------

To install, you need to run a USRP N200 with the custom fpga image and firmware, which are supplied for revision 4 of the device: n200_r4_fpga.bin, and usrp_n200_fw.bin. Refer the the Ettus Research UHD manual for instructions on installing firmware on the device. The source code is located in tarball UHD-003_005_002-phase_comparator.tar.gz.

The custom firmware simply provides decimated ADC samples on two channels, bypassing all the digital downconversion and filtering that would otherwise occur on the device (and screw up the phase precision in this process). 

Software
--------

You need Digital RF, a hdf5 based timestamped RF write library developed at Haystack. This depends on libhdf5-dev

> cd digital_rf/trunk
> ./configure 
> make 
> sudo make install
> sudo python setup.py install

You need gnuradio 3.7.6, and you need the out-of-tree module gr-drf. gr-drf contains a double precision digital downconverter that allows converting 5 MHz and 5 MHz + 1 kHz into baseband using sufficient precision to allow phase tracking with sub picosecond accuracy. The code is in dddc_impl.cc and it is written in C++ to ensure sufficient performance. Dropped packets are padded with zeros to keep time aligned. 

> cd gr-drf
> mkdir build
> cd build
> cmake ../
> make
> sudo make install
> sudo ldconfig

There is a recording script that can be run as a service to record phase of 5 MHz and 5 MHz + 1 kHz. The output goes into /data/phasecal. This directory needs to be user writeable. The data increases with about 100 MB/hour with the 100 Hz resolution. Old hourly directories can be deleted from the system while everything is running. The first 10 seconds of data are garbage, because this is used to determine DC offset. 

Starting the recorder:
> ./pcal_rec.py

A restart of the recorder will disrupt the continuity of the phase calibration due to resetting of the numerical oscillator. 

For outputting the cable delay as picoseconds, there is a command line tool that has several options summarized in the help.

./pcal_get_delay.py -i 100 -n --help
Usage: pcal_get_delay.py [options]

Options:
  -h, --help            show this help message and exit
  -d DIR, --dir=DIR     Directory. (default /data/phasecal)
  -i INTEGRATE, --integrate=INTEGRATE
                        Integration factor (default 100)
  -b BASELINE_TIME, --baseline_time=BASELINE_TIME
                        Time to use as reference for cable delay, unix
                        seconds. Default, now - 5 minutes
  -0 T0, --t0=T0        Start time in unix seconds, default: now-5 minutes
  -1 T1, --t1=T1        End time in unix seconds, default: now
  -p, --plot            plot relative time delay
  -o, --overview_plot   plot sparse overview plot
  -a, --ascii_out       output delays in ascii
  -n, --latest          Latest recorded delay

This tool requires three basic inputs: 
1) what timestamp is used to determine the reference phase (and cable delay)
2) what is the starting time for requesting values
3) what is the ending time for requesting values

For example, to output in ascii form the cable delay between unix second 1448903561 and 1448903571, with 1 second integration (100 Hz sample rate, integrate 100 samples), with the cable delay at unix second 1448903561 used as reference. The timestamp refers to the leading edge of the 1 second averaging window:

./pcal_get_delay.py -i 100 -b 1448903561 -0  1448903561 -1 1448903571 -a
# pcal out
# reference delay at 1448903561.00 (unix seconds), ref: 82761.43 (ps)
# integration 1.00 (seconds)
# time (unix seconds), delay (ps)
1448903561.00 0.00
1448903562.00 -0.27
1448903563.00 1.88
1448903564.00 5.64
1448903565.00 5.73
1448903566.00 5.74
1448903567.00 4.18
1448903568.00 1.48
1448903569.00 1.38
1448903570.00 4.96

A quick way to plot the last data recorded is to use -n. It is still advisable to supply a reference time each call is compared with a delay measured at the same time. 

> ./pcal_get_delay.py -n -b 1448903561
t0 1448904194.000 t1 1448904195.000 delay 3.271 reference time 1448903561.00

Values can also be plotted:
./pcal_get_delay.py -i 100 -b 1448903561 -0  1448903561 -1 1448903671 -p

Finally, there is a mode to sparsely go over a large amount of data (hard coded to 300 points evenly spread between t0 and t1). If one wants to enable to sparse mode, use the -o flag. This would e.g., plot 24 hours of data using 300 measurement points and 10 s resolution, plotting the result:

> ./pcal_get_delay.py -i 1000 -b 1448903561 -0  1448817161 -1 1448903621 -p -o

Warning: by default, the script will use now-60 seconds to determing reference time delay. This is not what you want in any operational measurement. 

Other
-----

The system is sensitive to temperature changes in cables and purity of the clock driving the system. When things are working correctly, ~ 1 picosecond stability is achieved. 
