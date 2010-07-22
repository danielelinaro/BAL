
AC_PREREQ([2.61])

# init autoconf, automake and libtool
AC_INIT([bal], [0.9.2], [daniele.linaro@unige.it])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_SRCDIR([src/balObject.h])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
LT_INIT

# operating system's information
AC_CANONICAL_HOST

# required programs
AC_PROG_AWK
AC_PROG_CXX
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# required header files
AC_CHECK_HEADERS([unistd.h])

# C language stuff
AC_C_CONST
AC_C_INLINE

# configuration macros
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_FILES([Makefile
		src/Makefile
		examples/Makefile
		python/setup.py
		python/pybal.h])

# required functions and types
AC_CHECK_FUNCS([floor pow sqrt])
AC_HEADER_STDBOOL
AC_TYPE_SIGNAL
AC_TYPE_SIZE_T
# required libraries
AC_SEARCH_LIBS([CVodeCreate], [sundials_cvode], [], [
echo "Can't find CVode library."
exit
])
AC_CHECK_LIB([sundials_cvode], [CVodeInit], [
CVODEVERSION=-DCVODE26
echo "CVode is version 2.6 or latter."], [])
AC_CHECK_LIB([sundials_cvode], [CVodeMalloc], [
CVODEVERSION=-DCVODE25
echo "CVode is version 2.5 or previous."], [])
AC_CHECK_LIB([sundials_nvecserial], [N_VNew_Serial], [], [
echo "Can't find NVecSerial library."
exit
])
AC_CHECK_LIB([hdf5], [H5Fcreate], [], [
echo "Can't find HDF5 library."
exit
])
AC_CHECK_LIB([hdf5_hl], [H5LTmake_dataset_double], [], [
echo "Can't find HDF5 library."
exit
])
AC_CHECK_LIB([boost_thread], [pthread_create], [], [
echo "Can't find boost thread library."
exit
])

# defines
LIB_CURRENT=0
LIB_REVISION=0
LIB_AGE=0
LIB_VER=$LIB_CURRENT:$LIB_REVISION:$LIB_AGE
LIB_TAG=$LIB_CURRENT.$LIB_REVISION.$LIB_AGE
AC_SUBST([LIB_VER])
AC_SUBST([LIB_TAG])
AC_SUBST([CVODEVERSION])


# conditional building
AC_ARG_WITH([examples],
	    [AS_HELP_STRING([--with-examples],[Build example programs])],
 	    [with_examples=yes],[echo "WARNING: not building the examples."])
AM_CONDITIONAL([EXAMPLES],[test "$with_examples" == "yes"])

echo "$host_os"
case "$host_os" in
  darwin*)
		EXTENSION="dylib"
    ;;
  linux*)
		EXTENSION="so"
    ;;
esac
AC_SUBST([EXTENSION])
		
AC_OUTPUT

