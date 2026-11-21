AC_PREREQ([2.70])
AC_INIT([biomc2], [2.0.1], [leomrtns@gmail.com])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([src/config.h])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_USE_SYSTEM_EXTENSIONS
AC_PROG_CC
AC_PROG_INSTALL
AC_SEARCH_LIBS([log], [m], [], [AC_MSG_ERROR([the maths library is required])])
AC_ARG_ENABLE([static-binary], [AS_HELP_STRING([--enable-static-binary], [link statically])])
AS_IF([test "x$enable_static_binary" = xyes], [LDFLAGS="$LDFLAGS -static"])
AC_ARG_ENABLE([optimization], [AS_HELP_STRING([--enable-optimization], [use optimised defaults unless CFLAGS was supplied])])
AS_IF([test "x$enable_optimization" = xyes && test -z "$ac_env_CFLAGS_set"], [CFLAGS="-O3"])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
