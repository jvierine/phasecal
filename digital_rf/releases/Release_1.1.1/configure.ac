#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#
#  Digital RF configure.ac file
#
#  Written by Bill Rideout (brideout@haystack.mit.edu)
#
#  $Id: configure.ac 522 2014-08-04 18:10:14Z brideout $

AC_PREREQ(2.59)
AC_INIT([Digital_rf], [1.1.1], [brideout@haystack.mit.edu])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for Hdf5 library

AC_CHECK_LIB([hdf5], [H5T_get_sign],,
  [AC_MSG_ERROR([Hdf5 C library must be installed first.  See http://www.hdfgroup.org/HDF5/])])
 

# distribute additional compiler and linker flags among Makefiles
# --> set and change these variables instead of CXXFLAGS or LDFLAGS (for user only)
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])

# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([source/Makefile])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
