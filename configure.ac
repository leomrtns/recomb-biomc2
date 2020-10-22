dnl Autoconf init, autoheader output location, config and prefix directories
AC_INIT([biomc2], [1.11], [leomrtns@gmail.com])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([src/config.h])
AC_PREFIX_DEFAULT(${HOME})

dnl automake initialization
AM_INIT_AUTOMAKE([-Wall -Werror])

dnl Basic compiler and related tools
AC_LANG_C
AC_PROG_CC
AC_PROG_INSTALL
AC_OPENMP

dnl necessary when using per-target flags (eg. foo_LDFLAGS); creates a "compile" shell script;
dnl when compiling will generate foo-file1.o ... instead of file1.o
AM_PROG_CC_C_O

dnl libraries for glib, gtk, etc (pkg-config --cflags --libs) 
dnl PKG_CHECK_MODULES(DEPS, gtk+-2.0 >= 2.2 glib-2.0 >= 2.2 gthread-2.0 cairo)
dnl PKG_CHECK_MODULES(DEPS, gtk+-2.0 >= 2.2 glib-2.0 >= 2.2 gthread-2.0)
PKG_CHECK_MODULES(DEPS, cairo)
dnl PKG_CHECK_MODULES([cairo], [cairo ])
AC_SUBST(DEPS_CFLAGS)
AC_SUBST(DEPS_LIBS)

dnl check for the math library
AC_CHECK_LIB([m],[main])

dnl Check for the 15 Standard C (ANSI-C/ISO-C) headers
AC_HEADER_STDC

dnl Check for additional headers 
AC_CHECK_HEADERS( \
    unistd.h \
    sys/time.h \
    sys/times.h \
)

AC_MSG_RESULT([....  configuration options specific to biomc2 ....])

AC_MSG_CHECKING([whether you want static binaries])
AC_ARG_ENABLE(static-binary,
    [  --enable-static-binary  static binaries that run on same arch without the libraries [[default=no]]],
    [ statbin_use="yes" ],
    [ statbin_use=""    ]
)
if test -n "${statbin_use}"; then
    AC_MSG_RESULT([yes])
		LDFLAGS="-static ${LDFLAGS}";
else
    AC_MSG_RESULT([no])
fi

AC_MSG_CHECKING([whether you want optimized binaries])
AC_ARG_ENABLE(optimization,
    [  --enable-optimization   optimization flags (may result in faster code) [[default=no]]],
    [ optim_use="yes" ],
    [ optim_use=""    ]
)
if test -n "${optim_use}"; then
    AC_MSG_RESULT([yes])
		CFLAGS="${CFLAGS} -std=gnu11 -funroll-loops -fomit-frame-pointer -finline-functions -O4"; # gnu89 retired
else
    AC_MSG_RESULT([no])
		CFLAGS="${CFLAGS} -std=gnu11 -g -W -Wall -Werror -Wno-uninitialized";
fi

AC_MSG_RESULT([....  end of configuration options specific to biomc2 ....])

dnl generate makefiles
AC_OUTPUT([ Makefile src/Makefile])

echo \
"---- 
Configuration parameters for biomc2:

 Source code location:  ${srcdir}
 Compiler:              ${CC}
 Compiler flags:        ${CFLAGS}
 Linker flags:          ${LDFLAGS}
 openMP multithreading: ${OPENMP_CFLAGS}
 Install path:          ${prefix}

----"

