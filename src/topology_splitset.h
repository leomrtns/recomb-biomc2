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

/*! \file topology_splitset.h
 *  \brief Functions that use only the split bipartitions of topologies -- treating them as unrooted usually.
 *
 *  Use a "splitset" structure that copies the bipartition information of all nodes (so that original trees are
 *  untouched) and then modifies this splitset -- which is in fact composed of two vectors of bipartitions like in the
 *  recombination program biomc2.
 */

#ifndef _biomcmc_topology_splitset_h_
#define _biomcmc_topology_splitset_h_

#include "bipartition.h"

typedef struct hungarian_struct* hungarian;
struct hungarian_struct
{
  int **cost, *col_mate; /*! \brief cost matrix, and col_mate[row] with column match for row */
  int size,  /*! \brief assignment size. Cost is a square matrix, so size should be an overestimate where "missing" nodes are added w/ cost zero */
      initial_cost, /*! \brief sum of lowest input cost values for each column. The hungarian method rescales them so that minimum per column is zero */
      final_cost;   /*! \brief our final cost is on rescaled cost matrix, therefore to restore the "classical" optimal cost one should sum it with initial_cost */
  int *unchosen_row, *row_dec, *slack_row, *row_mate, *parent_row, *col_inc, *slack; /* aux vectors */
};



typedef struct splitset_struct* splitset;

struct splitset_struct
{
  int size, spr, spr_extra, rf, hdist, hdist_reduced; /*! \brief spr, extra prunes for spr, rf distances and hdist=assignment cost */
  int n_g, n_s, n_agree, n_disagree;
  bipartition *g_split, *s_split, *agree, *disagree;
  int *leaf_idx, *leaf_recomb; /* original taxon -> active bit, and optional prune record */
  int cutoff;
  bipartition prune;
  hungarian h; /* hungarian method for solving the assignment between edges */
  bool match;  /*! \brief do we want to calculate the minimum cost assignment */
};

hungarian new_hungarian (int size);
void hungarian_reset (hungarian p);
void hungarian_update_cost (hungarian p, int row, int col, int cost);
void del_hungarian (hungarian p);
void hungarian_solve (hungarian p, int this_size);


splitset new_splitset_dSPR (int nleaves);
void del_splitset (splitset split);
int compare_splitset_bipartition_increasing (const void *a, const void *b);
int dSPR_topology_lowlevel (splitset split);
int rf_hdist_topology_lowlevel (splitset split, bool exit_at_rf);
#endif
