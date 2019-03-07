/*! 
 * \mainpage 
 * 
 *\section legal Legal Notice
 * 
 * Copyright (C) 2006	Leonardo de Oliveira Martins
 * 
 * leo at lbm ab a u-tokyo ac jp 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the 
 * Free Software Foundation, Inc., 51 Franklin Street, 
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*! \defgroup biomc2_global Global variables/functions
 */

/*! \file biomc2.h 
 *  \ingroup biomc2_global
 *  \brief Header file for biomc2.c 
 */

#ifndef _biomc2_h_
#define _biomc2_h_

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "config.h"

/* [ANSI C C89] */
#include <stdio.h> 
/* random number, searching, sorting [ANSI C C89] */
#include <stdlib.h>
/* string manipulation [ANSI C C89] */
#include <string.h>
/* Access to va_arg, va_list [ANSI C C89] */
#include <stdarg.h>
/* char operation functions (e.g. isspace() ), case convertion [ANSI C C89] */
#include <ctype.h>
/* standard math functions (e.g. exp() ) [ANSI C C89] */
#include <math.h>
/* speed profiling in CLOCKS_PER_SEC (e.g. clock() ) [ANSI C C89] */
#include <time.h>
/* system values checking at runtime (e.g. sysconf() ) [POSIX C] */
#include <unistd.h>
/* random seed (e.g. gettimeofday(), struct timeval) [POSIX C] */
#include <sys/time.h>
/* speed profiling in clock ticks (e.g. times() ) [POSIX.1 (or GNU extension?)] */ 
#include <sys/times.h>

/* Global constants */

#define EXP_1 2.718281828459045091

/*! \brief Program name and version (from autoconf) */
#define ProgramVersion PACKAGE_STRING 

#define true  1 /*!< Boolean TRUE  */
#define false 0 /*!< Boolean FALSE */

#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#define MOD(a)   (((a)>0)   ? (a) :(-a))

/*! \brief Maximum number of taxa the current implementation can handle 
	(owing to the split size limitation of 64 bits) */
#define MaxLeaves 64

/* Global Structures */

/*! \brief Mnemonic for boolean (char is smaller than int) */
typedef unsigned char bool;

/*! \brief Memory-safe malloc() function.
 *
 *  Allocates size bytes and returns a pointer to the allocated memory.
 *  An error message is thrown in case of failure.
 *  \param[in] size allocated size, in bytes
 *  \return pointer to newly allocated memory
 */
void *biomc2_malloc (size_t size);

/*! \brief Memory-safe realloc() function.
 *
 * Changes the size of the memory block pointed to by ptr to size bytes.
 * An error message is thrown in case of failure.
 * \param[in] size allocated size, in bytes
 * \param[in,out] ptr pointer to previously allocated memory
 * \return pointer to newly allocated memory
 */
void *biomc2_realloc (void *ptr, size_t size);

/*! \brief Memory-safe fopen() function.
 *
 * Opens the file whose name is the string pointed to by path and associates a
 * stream with it.
 * An error message is thrown in case of failure.
 * \param[in] path file name 
 * \param[in] mode opening mode ("r" for reading, "w" for writing, etc)
 * \result pointer to file stream
 */
FILE *biomc2_fopen (const char *path, const char *mode);

/*! \brief Prints error message and quits program.
 *
 * similar to fprintf (stderr, ...), but exits after printing the message
 * \param[in] template va_list following same format as printf()
 * \result exits program
 */
void biomc2_error (const char *template, ...);

/*! \brief Comparison between two integers.
 *
 * Comparison function used by sort(). Results in crescent order (from smaller
 * to larger).
 * \param[in] a pointer to integer
 * \param[in] b pointer to integer
 * \return 
 * - 1 if a>b
 * - -1 if a<b
 * - 0 if they are the same
 */
int compare_int (const void *a, const void *b);

/*! \brief Comparison between two doubles.
 *
 * Comparison function used by sort(). Results in crescent order (from smaller
 * to larger).
 * \param[in] a pointer to double
 * \param[in] b pointer to double
 * \return 
 * - 1 if a>b
 * - -1 if a<b
 * - 0 if they are the same
 */
int compare_double (const void *a, const void *b);

/*! \brief read file line-by-line (like homonymous function from GNU C library)
 *
 * This implementation is originally from the CvsGui project
 * (http://www.wincvs.org/). The explanation from the original file adapted to
 * our system  follows:
 * \verbatim 
   Read up to (and including) a newline ("\n") from STREAM into *LINEPTR
   and null-terminate it. *LINEPTR is a pointer returned from
   malloc (or NULL), pointing to *N characters of space.  It is realloc'd
   as necessary.  Return the number of characters read (not including the
   null terminator), or -1 on error or EOF.  
   \endverbatim
 */
int biomc2_getline (char **lineptr, size_t *n, FILE *stream);

#endif
