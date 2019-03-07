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

/*! \file nexus.c
 *  \ingroup nexus 
 *  \brief File handling functions for nexus format in general.
 */
#include "nexus.h"

/*! \brief Count balance of nexus comments (number of '['s minus number of ']'s). */
int   count_nested_nexus_comments (char *string);
/*! \brief Remove comments that start on this line (we need to guarantee that
 * they also finish on this line). */
char* remove_oneline_nexus_comments (char *string);

char*
remove_nexus_comments (char **string, size_t *stringsize, FILE *stream)
{
  int nested_comment = count_nested_nexus_comments (*string);
  while (nested_comment > 0) {
    if (biomc2_getline (string, stringsize, stream) == -1)  
      biomc2_error ( "Premature end of file when removing nexus comments\n");
    nested_comment += count_nested_nexus_comments (*string);
  }
  *string = remove_oneline_nexus_comments (*string);
  return (*string);
}

int 
count_nested_nexus_comments (char *string)
{
  char *s, *last=string+strlen (string);
  int count = 0;
  for (s = string; s <= last; s++) {
    if (*s == '[') count++;
    else if (*s == ']') count--;
  }
 return count;
} 

char* 
remove_oneline_nexus_comments (char *string) 
{
  char *s, *first, *last = string+strlen (string);
  int count = 1; /* we only use it if we find a "[" */
  for (s = string; (s <= last) && (*s != '['); s++);
  if (s>last) return (string);
  first = s++;
  while ((count > 0) && (s<=last)) {
    if ( *s == '[') count++;
    else if (*s == ']') count--;
    s++;
  }
	if (s<last) {
		*last = '\0';
		memmove (first, s, last - s + 1);
	}
  else return NULL;
  return string;
}

char*
lowercase_string (char *string)
{
  char *s, *last = string+strlen (string); 
  for (s = string; s <= last; s++)
    if (isupper (*s)) *s = tolower (*s);
  return string;
}

char*
uppercase_string (char *string)
{
  char *s, *last = string+strlen (string); 
  for (s = string; s <= last; s++)
    if (islower (*s)) *s = toupper (*s);
  return string;
}


char*
remove_space_from_string (char *string)
{
  char *last = strchr (string, '\0'); 
  char *s;

	for (s = last-1; s >= string; s--) {
		if (isspace (*s) || (*s == '>')) { memmove (s, s + 1, last - s); last--; }
	}
	return string;
}

