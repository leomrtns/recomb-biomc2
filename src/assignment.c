/*
 * This file is part of guenomu, a hierarchical Bayesian procedure to estimate the distribution of species trees based
 * on multi-gene families data.
 * Copyright (C) 2009  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * Guenomu is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file lowlevel.c
 *  \brief Lowest level basic functions, that should be available to all other modules.
 */

#include "topology_splitset.h"
/* The hungarian method below is copied from http://www.informatik.uni-freiburg.de/~stachnis/misc.html
 * The (edited) original message follows:
 *
 ** libhungarian by Cyrill Stachniss, 2004  Solving the Minimum Assignment Problem using the
 ** Hungarian Method.         ** This file may be freely copied and distributed! **
 **
 ** Parts of the used code was originally provided by the "Stanford GraphGase", but I made changes to this code.
 ** As asked by  the copyright node of the "Stanford GraphGase", I hereby proclaim that this file are *NOT* part of the
 ** "Stanford GraphGase" distrubition! */

void
hungarian_reset (hungarian p)
{
  int i, j;

  for (i = 0; i < p->size; i++) {
    p->col_mate[i] = p->unchosen_row[i] = p->row_dec[i] = p->slack_row[i] = p->row_mate[i] = p->parent_row[i] = p->col_inc[i] = p->slack[i] = 0;
    for (j = 0; j < p->size; j++) p->cost[i][j] = 0;
  }
  p->final_cost = 0;
}

hungarian
new_hungarian (int size)
{
  int i;
  hungarian p;

  p = (hungarian) biomc2_malloc (sizeof (struct hungarian_struct));
  p->size = size; /* n_rows = n_columns; if it's not, fill with zeroes (no cost) */
  p->cost = (int**) biomc2_malloc (size * sizeof (int*));
  for (i = 0; i < p->size; i++)
    p->cost[i] = (int*) biomc2_malloc (size * sizeof (int));
  /* edges would be assignment_matrix[ i * ncols + col_mate[i] ] = true; and other elems "false" (but we don't use the matrix notation) */
  p->col_mate     = (int*) biomc2_malloc (size * sizeof (int)); /* for a given row node, col_mate[row] is the assigned col node */
  p->unchosen_row = (int*) biomc2_malloc (size * sizeof (int));
  p->row_dec      = (int*) biomc2_malloc (size * sizeof (int));
  p->slack_row    = (int*) biomc2_malloc (size * sizeof (int));
  p->row_mate     = (int*) biomc2_malloc (size * sizeof (int));
  p->parent_row   = (int*) biomc2_malloc (size * sizeof (int));
  p->col_inc      = (int*) biomc2_malloc (size * sizeof (int));
  p->slack        = (int*) biomc2_malloc (size * sizeof (int));

  hungarian_reset (p);
  return p;
}

void
hungarian_update_cost (hungarian p, int row, int col, int cost)
{
  if (row >= p->size) return;
  if (col >= p->size) return;
  p->cost[row][col] = cost;
}

void
del_hungarian (hungarian p)
{
  int i;
  if (!p) return;
  if (p->cost) {
    for (i = p->size - 1; i >= 0; i--) if (p->cost[i]) free (p->cost[i]);
    free (p->cost);
  }
  free (p->col_mate); /* this is the important one, with i assigned to col_mate[i] */
  free (p->slack);
  free (p->col_inc);
  free (p->parent_row);
  free (p->row_mate);
  free (p->slack_row);
  free (p->row_dec);
  free (p->unchosen_row);
  free (p);
}

void
hungarian_solve (hungarian p, int this_size)
{
  int i, j, nrows = this_size, ncols = this_size, k, l, s, t, q, unmatched;
  p->final_cost = p->initial_cost = 0;

  if (this_size > p->size) { p->final_cost = -1; return; } /* we don't call biomcmc_error(), but it *is* an error! */

  for (l = 0; l < ncols; l++) { // Begin subtract column minima in order to start with lots of zeroes 12
    s = p->cost[0][l];
    for (k = 1; k < nrows; k++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->initial_cost += s; /* this should be added to final_cost to have classical assignment cost; here we distinguish them */
    if (s!=0)	for (k = 0; k < nrows; k++) p->cost[k][l] -= s;
  } // End subtract column minima in order to start with lots of zeroes 12

  /*for (i=0;i<nrows; i++) {
    for (j=0;j<ncols; j++) printf ("%4d ", p->cost[i][j]);
    printf (" DEBUGcost\n");
  }  // DEBUG */

  // Begin initial state 16
  t=0;
  for (l = 0; l < ncols; l++)  { // n => num_cols
    p->row_mate[l]= -1;
    p->parent_row[l]= -1;
    p->col_inc[l]=0;
    p->slack[l]= 0x7FFFFFFF;
  }
  for (k = 0; k < nrows; k++) { // m => num_rows
    s = p->cost[k][0];
    for (l = 1; l < ncols; l++) if (p->cost[k][l] < s) s = p->cost[k][l];
    p->row_dec[k]=s;
    for (l = 0; l < ncols; l++) if ((s==p->cost[k][l]) && (p->row_mate[l] < 0)) {
      p->col_mate[k] = l;
      p->row_mate[l] = k;  // fprintf(stderr, "matching col %d==row %d\n",l,k);
      goto row_done;
    }
    p->col_mate[k] = -1;  // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    p->unchosen_row[t++] = k;
row_done:
    ;
  }
  // End initial state 16

  // Begin Hungarian algorithm 18
  if (t==0)    goto done;
  unmatched=t;
  while (1) {
    q=0; // fprintf(stderr, "Matched %d rows.\n",m-t);
    while (1)	{
      while (q<t) {
         { // Begin explore node q of the forest 19
          k = p->unchosen_row[q];
          s=p->row_dec[k];
          for (l=0;l<ncols;l++) if (p->slack[l]) {
            int del;
            del = p->cost[k][l] - s + p->col_inc[l];
            if (del < p->slack[l]) {
              if (del==0) {
                if (p->row_mate[l]<0)  goto breakthru;
                p->slack[l]=0;
                p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n", t,row_mate[l],l,k);
                p->unchosen_row[t++]=p->row_mate[l];
              }
              else { p->slack[l]=del; p->slack_row[l]=k; }
            }
          }
         } // End explore node q of the forest 19
        q++;
      }

      // Begin introduce a new zero into the matrix 21
      s = 0x7FFFFFFF;
      for (l = 0;l < ncols; l++) if (p->slack[l] && p->slack[l] < s) s = p->slack[l];
      for (q = 0; q < t; q++) p->row_dec[ p->unchosen_row[q] ] += s;
      for (l = 0; l < ncols; l++) if (p->slack[l]) {
        p->slack[l]-=s;
        if (p->slack[l]==0) {  // Begin look at a new zero 22
          k = p->slack_row[l]; // fprintf(stderr, "Decreasing uncovered elements by %d produces zero at [%d,%d]\n", s,k,l);
          if (p->row_mate[l]<0)  {
            for (j=l+1;j<ncols;j++)  if (p->slack[j]==0) p->col_inc[j]+=s;
            goto breakthru;
          }
          else {
            p->parent_row[l]=k; // fprintf(stderr, "node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
            p->unchosen_row[t++]=p->row_mate[l];
          }
        } // End look at a new zero 22

      }
      else  p->col_inc[l]+=s;
      // End introduce a new zero into the matrix 21
    }
breakthru:
    // fprintf(stderr, "Breakthrough at node %d of %d!\n",q,t);
    while (1)	{    // Begin update the matching 20
      j=p->col_mate[k];
      p->col_mate[k]=l;
      p->row_mate[l]=k; // fprintf(stderr, "rematching col %d==row %d\n",l,k);
      if (j<0)   break;
      k=p->parent_row[j];
      l=j;
    }    // End update the matching 20
    if (--unmatched==0)	goto done;
    // Begin get ready for another stage 17
    t=0;
    for (l=0;l<ncols;l++) {
      p->parent_row[l]= -1;
      p->slack[l]=0x7FFFFFFF;
    }
    for (k=0;k<nrows;k++) if (p->col_mate[k]<0) p->unchosen_row[t++]=k; // fprintf(stderr, "node %d: unmatched row %d\n",t,k);
    // End get ready for another stage 17
  }
done:
  // Begin doublecheck the solution 23
  for (k = 0; k < nrows; k++) for (l=0;l<ncols;l++) if (p->cost[k][l] < p->row_dec[k] - p->col_inc[l]) { p->final_cost = -1; printf ("\n**\n"); return;}
  for (k = 0; k < nrows; k++) {
    l=p->col_mate[k];
    if ((l < 0) || (p->cost[k][l] != p->row_dec[k] - p->col_inc[l])) { p->final_cost = -1; return; }
  }
  k=0;
  for (l=0;l<ncols;l++) if (p->col_inc[l])  k++;
  if (k>nrows) { p->final_cost = -1; return; }
  // End doublecheck the solution 23
  // End Hungarian algorithm 18

  for (k = 0; k < nrows; ++k) for (l = 0; l < ncols; ++l) p->cost[k][l] = p->cost[k][l] - p->row_dec[k] + p->col_inc[l];
  for (i = 0; i < nrows; i++) p->final_cost += p->row_dec[i];
  for (i = 0; i < ncols; i++) p->final_cost -= p->col_inc[i]; // fprintf(stderr, "Cost is %d\n",cost);
}
