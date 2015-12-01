# simple example Makefile for rf_write_hdf5 C code
#
# $Id: Makefile 386 2014-05-23 14:44:59Z brideout $
all: test_rf_write_hdf5 benchmark_rf_write_hdf5 libdigital_rf.a

CFLAGS=-fPIC

test_rf_write_hdf5: test_rf_write_hdf5.o rf_write_hdf5.o
	gcc -Wall $(CFLAGS) -I. -I/usr/local/include -o test_rf_write_hdf5 -g rf_write_hdf5.c test_rf_write_hdf5.c -lhdf5 -lm

libdigital_rf.a: rf_write_hdf5.o 
	gcc -Wall $(CFLAGS) -I. -I/usr/local/include -c rf_write_hdf5.c -lhdf5 -lm
	rm -rf $@ 
	ar cq $@ rf_write_hdf5.o 

benchmark_rf_write_hdf5: benchmark_rf_write_hdf5.o rf_write_hdf5.o
	gcc -Wall $(CFLAGS) -I. -I/usr/local/include -o benchmark_rf_write_hdf5 -g rf_write_hdf5.c benchmark_rf_write_hdf5.c -lhdf5 -lm

install: libdigital_rf.a
	cp libdigital_rf.a /usr/lib/
	cp digital_rf.h /usr/include/digital_rf.h