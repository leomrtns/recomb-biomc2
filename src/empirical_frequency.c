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

/*! \file empirical_frequency.c  
 *  \brief histogram of vectors, ordered by frequency. Also calculates MAP (modal) values.
 *
 * Sorts a vector of integers by their frequencies, preserving their original indexes. It is a simple extension to qsort
 * where the original order can be reconstructed, or still a key/value sorting.
 */

#include "empirical_frequency.h"

int compare_empfreq_element_increasing (const void *a, const void *b);
int compare_empfreq_element_decreasing (const void *a, const void *b);

int
compare_empfreq_element_decreasing (const void *a, const void *b)
{
  int result;
  empfreq_element *e_a = (empfreq_element *) a;
  empfreq_element *e_b = (empfreq_element *) b;
  result =  e_b->freq - e_a->freq;
  if (result) return result;
  /* break ties by position */
  return (e_a->idx - e_b->idx);
}

int
compare_empfreq_element_increasing (const void *a, const void *b)
{
  empfreq_element *e_a = (empfreq_element *) a;
  empfreq_element *e_b = (empfreq_element *) b;
  return (e_a->freq - e_b->freq);
}

void
sort_empfreq_decreasing (empfreq ef)
{
  qsort (ef->i, ef->n, sizeof (empfreq_element), compare_empfreq_element_decreasing);
}

void
sort_empfreq_increasing (empfreq ef)
{
  qsort (ef->i, ef->n, sizeof (empfreq_element), compare_empfreq_element_increasing);
}

empfreq
new_empfreq (int n_elements)
{
  empfreq ef;
  int i;

  ef = 
  (empfreq) biomc2_malloc (sizeof (struct empfreq_struct));
  ef->n = n_elements;
  ef->min = 0;
  ef->max = n_elements - 1;
  ef->i =
  (empfreq_element*) biomc2_malloc (n_elements * sizeof (empfreq_element));

  for (i=0; i< n_elements; i++) {
    ef->i[i].freq = 0;
    ef->i[i].idx  = i;
  }

  return ef;
}

void
del_empfreq (empfreq ef)
{
  if (ef) {
    if (ef->i) free (ef->i);
    free (ef);
  }
}

empfreq
new_empfreq_from_int_weighted (int *vec, int n, int *weight)
{
  int i, distinct_values = 1, new_n = 0, nonzero[n];
  empfreq e_idx, e_count;

  /* excludes elements of vec with weight=0 */
  for (i=0; i < n; i++) if (weight[i]) { nonzero[new_n++] = i; }


  if (!new_n) biomc2_error ("vector of empirical frequencies with all freqs=0");

  e_idx = new_empfreq (new_n);
  for (i=0; i < new_n; i++) {
    /* trick to make it sort by value and not by frequency (the "real" sort comes below) */
    e_idx->i[i].idx  = weight[ nonzero[i] ];
    e_idx->i[i].freq =    vec[ nonzero[i] ];
  }
  sort_empfreq_increasing (e_idx); /* equiv. to qsort (vec, n, sizeof (int), compare_int) but preserving weights. */

  for (i=1; i < new_n; i++) if (e_idx->i[i].freq != e_idx->i[i-1].freq) distinct_values++;

  e_count = new_empfreq (distinct_values);

  for (distinct_values = 0, i = 0; i < new_n-1; i++) {
    e_count->i[distinct_values].freq += e_idx->i[i].idx;
    e_count->i[distinct_values].idx = e_idx->i[i].freq;
    if (e_idx->i[i].freq != e_idx->i[i+1].freq) distinct_values++;
  }
  e_count->i[distinct_values].freq += e_idx->i[i].idx;
  e_count->i[distinct_values].idx = e_idx->i[i].freq;

  e_count->min = e_count->i[0].idx;
  e_count->max = e_count->i[distinct_values].idx;

  sort_empfreq_decreasing (e_count); /* Now it will sort from largest to smallest freq (sum of weights) */
  del_empfreq (e_idx);
  return e_count;
}

empfreq
new_empfreq_from_int (int *vec, int n)
{
  int i, weight[n];
  for (i=0; i<n; i++) weight[i]=1;
  return new_empfreq_from_int_weighted (vec, n, weight);
}

int
find_mode_int_weighted (int *vec, int n, int *weight)
{
  empfreq e_count;
  int map;
  e_count = new_empfreq_from_int_weighted (vec, n, weight);
  map = e_count->i[0].idx;
  del_empfreq (e_count);
  return map;
}

int
find_mode_int (int *vec, int n)
{
  int i, weight[n];
  for (i=0; i<n; i++) weight[i]=1;
  return find_mode_int_weighted (vec, n, weight);
}

