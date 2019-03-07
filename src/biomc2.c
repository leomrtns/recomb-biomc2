/* Copyright (C) 2006  Leonardo de Oliveira Martins
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

/*! \file biomc2.c 
 *  \ingroup biomc2_global 
 *  \brief Global functions
 *
 *  This file contains basic functions that should be available to all other
 *  modules. 
 */

#include "biomc2.h"

/* error-safe memory allocation functions */

void *
biomc2_malloc (size_t size)
{
  void *value = malloc (size);
  if (value == NULL)
    biomc2_error ( "biomc2_malloc ");
  return value;
}

void *
biomc2_realloc (void *ptr, size_t size)
{
  void *value = (void *) realloc ((void *)ptr, size);
  if (value == NULL)
    biomc2_error ( "biomc2_realloc ");

//  fprintf(stderr, "realloc: 0x%08X => 0x%08X [%d]\n", ptr, value, size);
  return value;
}

FILE *
biomc2_fopen (const char *path, const char *mode)
{
	FILE *fp = fopen (path, mode);
	if (fp == NULL) 
		biomc2_error ( "biomc2_fopen_read (file %s) ", path);
	return fp;
}

void
biomc2_error (const char *template, ...)
{
	va_list ap;

	fprintf (stderr, "biomc2 error: ");
	va_start (ap, template);
	vfprintf (stderr, template, ap);
	va_end (ap);
	fprintf (stderr, "\n");
	exit (EXIT_FAILURE);
}

/* increasing order */
int
compare_int (const void *a, const void *b)
{
	return (*(int *) a - *(int *) b);
}

int
compare_double (const void *a, const void *b)
{
	if ((double *) a > (double *) b) return 1;
	if ((double *) a < (double *) b) return -1;
	return 0;
}

/* \brief size, in bytes, when extending the buffer of biomc2_getline() */
#define MIN_CHUNK 64
int
biomc2_getline (char **lineptr, size_t *n, FILE *stream)
{
	int nchars_avail;		/* Allocated but unused chars in *LINEPTR.  */
	char *read_pos;		  /* Where we're reading into *LINEPTR. */

	if (!lineptr) biomc2_error ("NULL pointer sent to biomc2_getline() as target string");
	if (!n)       biomc2_error ("string length unavailable to biomc2_getline()");
	if (!stream)  biomc2_error ("lack of input file in biomc2_getline()");

	if (!(*lineptr)) {
		*n = MIN_CHUNK;
		*lineptr = (char *) biomc2_malloc (*n);
	}
	nchars_avail = *n;
	read_pos = *lineptr;

	for (;;) {
		register int c = getc (stream);

		/* We always want at least one char left in the buffer, since we
		 * always (unless we get an error while reading the first char)
		 * NUL-terminate the line buffer.  
		 */
		if ((*lineptr + *n) != (read_pos + nchars_avail)) 
			biomc2_error ("problem setting string size in biomc2_getline()");
		if (nchars_avail < 2) {
			if (*n > MIN_CHUNK) (*n) *= 2;
			else (*n) += MIN_CHUNK;

			nchars_avail = *n + *lineptr - read_pos;
			*lineptr = (char *) biomc2_realloc (*lineptr, *n);
			read_pos = *n - nchars_avail + *lineptr;
			if ((*lineptr + *n) != (read_pos + nchars_avail)) 
				biomc2_error ("problem setting string size in biomc2_getline()");
		}

		if (ferror (stream)) return -1;

		if (c == EOF) {
			/* Return partial line, if any.  */
			if (read_pos == *lineptr) return -1;
			else break;
		}

		if (c == '\r') c = '\n';

		*read_pos++ = c;
		nchars_avail--;

		if (c == '\n') break;
	}

	/* Done - NUL terminate and return the number of chars read.  */
	*read_pos = '\0';
	return (read_pos - (*lineptr));
}

