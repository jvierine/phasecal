# $Id: Makefile.am 818 2015-09-16 16:19:41Z brideout $
#
# Makefile.am for Digital RF
#
# Written by Bill Rideout (brideout@haystack.mit.edu)

ACLOCAL_AMFLAGS = -I m4

include_HEADERS = digital_rf.h

lib_LTLIBRARIES = libdigital_rf.la

libdigital_rf_la_SOURCES = rf_write_hdf5.c

LDADD = -ldigital_rf

# test programs only
check_PROGRAMS = test_rf_write_hdf5 benchmark_rf_write_hdf5 example_rf_write_hdf5

