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

/*! \file nexus_alignment.c
 *  \ingroup nexus 
 *  \brief File handling functions for sequence data in nexus format.
 */

#include "nexus.h"

/* local function prototypes */
/*! \brief Allocates space for new alignment struct. */
nexus_alignment new_nexus_alignment (int ntax, int nchar);
/*! \brief Reads one line of alignment (sequence data for one taxa). */
void  read_interleaved_nexus_sequence (char *line, nexus_alignment align);
/*! \brief Allocates and stores taxa name. */ 
char* new_taxlabel (char *name);
/*! \brief Reads partition information from alignment. */
void read_assumptions_nexus_sequence (char *line, nexus_alignment align);

nexus_alignment
read_nexus_alignment_file (char *seqfilename)
{
  nexus_alignment align=NULL;
  FILE *seqfile;
  char *line = NULL, *line_read = NULL, *needle_tip;
	bool option_begin_data   = false, 
			 option_begin_matrix = false, 
			 option_end_matrix   = false, 
			 option_interleave   = false,
			 option_begin_assumptions  = false,
			 option_end_assumptions = false,
			 option_read_complete = false;
  size_t linelength = 0;
  int ntax=0, nchar=0;

  seqfile = biomc2_fopen (seqfilename, "r");
  biomc2_getline (&line_read, &linelength, seqfile);

  if (!strcasestr (line_read, "NEXUS")) biomc2_error ( "Not a Nexus sequence file\n");

	/* the variable *line should point always to the same value (no line++ or
	 * alike) */
  while (biomc2_getline (&line_read, &linelength, seqfile) != -1) {
    line = remove_nexus_comments (&line_read, &linelength, seqfile);
    /* do not parse empty lines ( = "\n\0") */	
    if (strlen (line) > 1) { 
    
      /* do not do anything untill 'BEGIN DATA' block */
      if ((!option_begin_data) && (strcasestr (line, "BEGIN DATA"))) option_begin_data = true;

      /* read ntax and nchar info, and wait for 'MATRIX' block */
      else if (!option_begin_matrix) {
        if (strcasestr (line, "MATRIX")) option_begin_matrix = true;
        else if (strcasestr (line, "INTERLEAVE")) option_interleave = true;
        /* needle_tip will point to first character of 'DIMENSIONS', removing
         * indenting spaces */
        else if ( (needle_tip = strcasestr (line, "DIMENSIONS")) ) {
          needle_tip += strlen("DIMENSIONS");
          needle_tip = uppercase_string (needle_tip);
          sscanf (needle_tip, " NTAX = %d NCHAR = %d ", &ntax, &nchar);
          align = new_nexus_alignment (ntax, nchar);
        } // else if (needle_tip)
      }

      /* entering 'MATRIX' block with sequence data */
      else if (!option_read_complete){
        /* semi-colon may be at the last line of seq data, or at new one */
        if ( (needle_tip = strstr (line, ";")) ) {
          option_read_complete = true;
          if (needle_tip > line) { 
            if (option_interleave) read_interleaved_nexus_sequence (line, align);
            //else read_sequencial_nexus_sequence (line, align);
            else biomc2_error ( "seq data should be interleaved\n");
          }
        } // if (needle_tip = strstr())
        else {	
          if (option_interleave) read_interleaved_nexus_sequence (line, align);
          //else read_sequencial_nexus_sequence (line, align);
          else biomc2_error ( "seq data should be interleaved\n");
        }
      }

      /* after finish reading matrix data, we need to wait for END keyword */
      else if ((!option_end_matrix) && (strcasestr (line, "END")))
        option_end_matrix = true;

      /* wait for assumptions block */
      else if ((!option_begin_assumptions) && (strcasestr (line, "BEGIN ASSUMPTION")))
        option_begin_assumptions = true; 

      /* once in assumptions block, just read CHARSET (gene segments) and order
       * segments after END*/ 
      else if (!option_end_assumptions) {
        if (strcasestr (line, "CHARSET") && (needle_tip = strcasestr (line, "="))) {
					needle_tip++;
					read_assumptions_nexus_sequence (needle_tip, align);
				}
				else if (strcasestr (line, "END")) {
					option_end_assumptions = true;
					qsort (align->charset_start, align->n_charset, sizeof (int), compare_int);
					qsort (align->charset_end, align->n_charset, sizeof (int), compare_int);
				}
			}
			
			/* other commands, at the end of file, come here (empty now) */
			
		} // if (line)
	} //while (biomc2_getline)
	fclose (seqfile);
	if (line) free (line);
	return align;
}

void
print_nexus_alignment (nexus_alignment align) /*globalfunction */
{
	int i, row, columns = 80;

	for (i = 0; i < align->ntax; i++)
	 {
		printf (">%s\n", align->taxlabel[i]);
		for (row = 0; row < align->nchar; row += columns)
			printf ("%.*s\n", columns, align->character[i] + row);
	 }
}

nexus_alignment
new_nexus_alignment (int ntax, int nchar) /* localfunction */
{
  nexus_alignment align;
  int i;

  align = 
  (nexus_alignment) biomc2_malloc (sizeof (struct nexus_alignment_struct));
  align->ntax = ntax;
  align->nchar = nchar;

  align->n_charset = 0;
  align->charset_start = NULL; /* (re)allocated by read_assumptions() */
  align->charset_end   = NULL; /* (re)allocated by read_assumptions() */


  align->taxlabel = (char **) biomc2_malloc (ntax * sizeof (char*));
  for (i=0; i < ntax; i++) align->taxlabel[i] = NULL;
  /* align->taxlabel[i] is allocated by new_taxlabel function to economize
   * bites ;) */ 

  align->character = (char **) biomc2_malloc (ntax * sizeof (char*));
  for (i=0; i < ntax; i++) {
    align->character[i] = (char *) biomc2_malloc ((nchar+1) * sizeof (char));
    align->character[i][0] = '\0'; 
  }

  align->taxlabel_hash = new_hashtable (ntax);

  return align;
}

void 
del_nexus_alignment (nexus_alignment align) /* globalfunction */
{
	int i;

  if (align) {
    if (align->character) { 
      for (i=align->ntax-1; i>=0; i--) 
        if (align->character[i]) free (align->character[i]);
      free (align->character);
    }
		/* align->taxlabel is referred by chain_data and hence cannot be deleted
		 * here. */
		/*
		if (align->taxlabel) { 
			for (i=align->ntax-1; i>=0; i--) 
				if (align->taxlabel[i])  free (align->taxlabel[i]);
			free (align->taxlabel);
		}
		 */
    del_hashtable (align->taxlabel_hash);
    if (align->charset_start) free (align->charset_start);
    if (align->charset_end) free (align->charset_end);
    free (align);
  }
}

void
read_interleaved_nexus_sequence (char *line, nexus_alignment align) /* localfunction */
{
  char seqname[MAX_NAME_LENGTH]="", *last = NULL;
	int namelength; 
	int position;
	namelength = strcspn (line, " \t\n");
	if (namelength < MAX_NAME_LENGTH)
		strncpy (seqname, line, namelength);
	else biomc2_error ( "taxa name too long\n");
	position = lookup_hashtable (align->taxlabel_hash, seqname);
	if (position < 0) {
		position = binary_search_next_available (align->taxlabel, align->ntax);
		align->taxlabel[position] = new_taxlabel (seqname);
		insert_hashtable (align->taxlabel_hash, seqname, position);
	}
	line = remove_space_from_string (line+namelength);
  if ((last = strchr (line, ';')))
    strncat (align->character[position], line, last - line - 1);
	else
		strncat (align->character[position], line, strlen (line));
}

char *
new_taxlabel (char *name)
{
	char *string;
	string = (char*) biomc2_malloc (sizeof (char) * (strlen (name)+1));
	strcpy (string, name);
	return string;
}

void
read_assumptions_nexus_sequence (char *line, nexus_alignment align)
{
	int start, end;

	align->n_charset++;
	align->charset_start = 
	(int*) biomc2_realloc ((int*) align->charset_start, sizeof (int) * align->n_charset);
	align->charset_end = 
	(int*) biomc2_realloc ((int*) align->charset_end, sizeof (int) * align->n_charset);
	
	if (sscanf (line, " %d - %d ;", &start, &end) != 2) 
		biomc2_error ( "problem in ASSUMPTIONS block\n");
  
//  printf ("entering assumpt with start = %d and end = %s\n", start, end);
	
	start--; end--;
	/* error handling */
	if (start < 0) biomc2_error ( "CHARSET start < 1\n");
	if (start >= align->nchar)  biomc2_error ( "CHARSET start > NCHAR\n");
	if (end <= start) biomc2_error ( "CHARSET end <= CHARSET start\n");
	if (end >= align->nchar)  biomc2_error ( "CHARSET end > NCHAR\n");
	
	align->charset_start[align->n_charset-1] = start;
	align->charset_end[align->n_charset-1] = end;
}

