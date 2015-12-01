# $Id: Makefile.am 794 2015-07-07 17:51:49Z brideout $
#
# Makefile.am for Digital RF
#
# Written by Bill Rideout (brideout@haystack.mit.edu)

ACLOCAL_AMFLAGS = -I m4

include_HEADERS = digital_rf.h

lib_LTLIBRARIES = libdigital_rf.la

libdigital_rf_la_SOURCES = rf_write_hdf5.c

LDADD = -ldigital_rf

docdir = $(prefix)/$(PACKAGE)/doc
dist_doc_DATA = documents/Digital_RF.pdf


# test programs only
check_PROGRAMS = test_rf_write_hdf5 benchmark_rf_write_hdf5

