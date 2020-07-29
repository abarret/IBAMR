## ---------------------------------------------------------------------
##
## Copyright (c) 2014 - 2014 by the IBAMR developers
## All rights reserved.
##
## This file is part of IBAMR.
##
## IBAMR is free software and is distributed under the 3-clause BSD
## license. The full text of the license can be found in the file
## COPYRIGHT at the top level directory of IBAMR.
##
## ---------------------------------------------------------------------

AC_DEFUN([CONFIGURE_SAMRAI],[
echo
echo "==================================="
echo "Configuring required package SAMRAI"
echo "==================================="

PACKAGE_SETUP_ENVIRONMENT

AM_CONDITIONAL([SAMRAI2D_ENABLED],true)
AM_CONDITIONAL([SAMRAI3D_ENABLED],true)
PACKAGE_RESTORE_ENVIRONMENT
])
