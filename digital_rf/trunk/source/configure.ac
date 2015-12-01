#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
#  Digital RF configure.ac file
#
#  Written by Bill Rideout (brideout@haystack.mit.edu)
#
#  $Id: configure.ac 787 2015-07-07 15:02:14Z brideout $

AC_PREREQ(2.59)
AC_INIT([Digital_rf], [1.1.3], [brideout@haystack.mit.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_PROG_AR
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
