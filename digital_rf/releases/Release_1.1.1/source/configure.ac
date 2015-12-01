#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
#  Digital RF configure.ac file
#
#  Written by Bill Rideout (brideout@haystack.mit.edu)
#
#  $Id: configure.ac 390 2014-05-23 15:21:00Z brideout $

AC_PREREQ(2.59)
AC_INIT([Digital_rf], [1.0], [brideout@haystack.mit.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for Hdf5 library

AC_CHECK_LIB([hdf5], [H5T_get_sign],,
  [AC_MSG_ERROR([Hdf5 C library must be installed first.  See http://www.hdfgroup.org/HDF5/])])
 


# Checks for library functions.

AC_CONFIG_FILES([
    Makefile
])
AC_OUTPUT
