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

#include "topology_splitset.h"

/*! \brief function used by qsort for a vector of bipartitions (from smaller to larger) */
int compare_splitset_bipartition_increasing (const void *a1, const void *a2);
/*! \brief function used by qsort for a vector of bipartitions (from larger to smaller) */
int compare_splitset_bipartition_decreasing (const void *a1, const void *a2);
/*! \brief copies bipartion information of internal nodes to splitset, "unroots" the bipartitions and order the resulting vectors */
int rf_hdist_topology_lowlevel (splitset split, bool exit_at_rf);
int dSPR_topology_lowlevel (splitset split);
/*! \brief vector with pointers to bipartitions that are identical on both trees */
void split_create_agreement_list (splitset split);
void split_remove_agree_edges (splitset split, bipartition *b, int *nb);
void split_remove_duplicates (bipartition *b, int *nb);
void split_compress_agreement (splitset split);
void split_create_disagreement_list (splitset split);
void split_disagreement_assign_match (splitset split);
void split_find_small_disagreement (splitset split);
void split_remove_small_disagreement (splitset split);
void split_minimize_subtrees (splitset split);
void split_remove_redundant_bit (splitset split, int id);
void split_replace_bit (splitset split, int to, int from);
void split_new_size (splitset split, int size, bool update_bipartitions);
void split_swap_position (bipartition *b, int i1, int i2);

splitset
new_splitset (int nleaves)
{
  splitset split;
  int i;

  split = (splitset) biomc2_malloc (sizeof (struct splitset_struct));
  split->size = nleaves - 1;
  split->disagree = split->agree = NULL; /* used only when calculating dSPR */
  split->n_agree = split->n_disagree = 0;
  split->leaf_idx = split->leaf_recomb = NULL;
  split->prune = NULL;
  split->h = NULL; /* hungarian method for bipartite matching; unused unless for tree distance */
  split->match = false;
  split->spr = split->spr_extra = split->rf = split->hdist = split->hdist_reduced = 0;

  split->g_split = (bipartition*) biomc2_malloc (split->size * sizeof (bipartition));
  split->s_split = (bipartition*) biomc2_malloc (split->size * sizeof (bipartition));
  split->g_split[0] = new_bipartition (nleaves);
  split->s_split[0] = new_bipartition (nleaves);
  for (i = 1; i < split->size; i++) {
    split->g_split[i] = new_bipartition_from_bipsize (split->g_split[0]->n); /* use same bipsize */
    split->s_split[i] = new_bipartition_from_bipsize (split->s_split[0]->n);
  }
  split->n_g = split->n_s = 0; /* number of elements actively in use */

  return split;
}

splitset
new_splitset_dSPR (int nleaves)
{
  int i;
  splitset split = new_splitset (nleaves);

  split->agree    = (bipartition*) biomc2_malloc (split->size * sizeof (bipartition));
  split->disagree = (bipartition*) biomc2_malloc (split->size * sizeof (bipartition));
  split->agree[0]    = new_bipartition (nleaves); // this bipsize will be recycled below
  split->disagree[0] = new_bipartition (nleaves);
  for (i = 1; i < split->size; i++)               split->agree[i]    = new_bipartition_from_bipsize (split->agree[0]->n);
  for (i = 1; i < split->size; i++) split->disagree[i] = new_bipartition_from_bipsize (split->disagree[0]->n);
  split->prune = new_bipartition_from_bipsize (split->disagree[0]->n);

  /* Allocate the assignment matrix only when a comparison has unmatched edges. */
  split->match = false; /* do we want to calculate the assignment matching cost (using hungarian() )? */

  return split;
}

void
del_splitset (splitset split)
{
  int i;
  if (!split) return;

  del_bipartition (split->prune);
  if (split->disagree) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->disagree[i]);
    free (split->disagree);
  }
  if (split->agree) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->agree[i]);
    free (split->agree);
  }
  if (split->g_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->g_split[i]);
    free (split->g_split);
  }
  if (split->s_split) {
    for (i = split->size - 1; i >= 0; i--) del_bipartition (split->s_split[i]);
    free (split->s_split);
  }
  del_hungarian (split->h);
  free (split);
}

int
compare_splitset_bipartition_increasing (const void *a1, const void *a2)
{ /* similar to bipartition_is_larger() */
  bipartition *b1 = (bipartition *) a1;
  bipartition *b2 = (bipartition *) a2;
  int i;

  if ((*b1)->n_ones > (*b2)->n_ones) return 1;
  if ((*b1)->n_ones < (*b2)->n_ones) return -1;

  for (i = (*b1)->n->ints - 1; (i >= 0) && ((*b1)->bs[i] == (*b2)->bs[i]); i--); /* find position of distinct bipartition elem*/
  if (i < 0) return 0; /* identical bipartitions */
  if ((*b1)->bs[i] > (*b2)->bs[i]) return 1;
  else return -1;
}

int
compare_splitset_bipartition_decreasing (const void *a1, const void *a2)
{
  return compare_splitset_bipartition_increasing (a2, a1);
}

int
rf_hdist_topology_lowlevel (splitset split, bool exit_at_rf)
{
  split->hdist_reduced = split->hdist = split->rf = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits);
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits);
  split_create_agreement_list (split);  // vector of identical bipartitions
  // Importantly, here we do NOT call split_compress_agreement()
  split->rf = split->n_g + split->n_s;
  if (exit_at_rf) return split->rf;
  if (!split->rf) return 0; // all edges were in agreement
  split->match = true; // only calculate hdist_reduced if match == True (first time)
  split_create_disagreement_list (split); // vector of smallest disagreements
  split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist_reduced after first time) */
  split->hdist = split->hdist_reduced;
  return split->hdist; /* do not calculate SPR, exit now */
}

int
dSPR_topology_lowlevel (splitset split)
{
  bool mismatch = true;
  split->match = true;
  split->hdist_reduced = split->spr = split->spr_extra = 0;
  split->n_agree = split->n_disagree = 0;
  bipsize_resize (split->disagree[0]->n, split->g_split[0]->n->bits);
  bipsize_resize (split->agree[0]->n,    split->g_split[0]->n->bits);

  while (mismatch) {
    split_create_agreement_list (split);  // vector of identical bipartitions
    split_compress_agreement (split);     // iterative replacement of cherry by new leaf

    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
    if (!mismatch) return split->spr;

    split_create_disagreement_list (split); // vector of smallest disagreements
    split_disagreement_assign_match (split); /* assignment matching between edges using hungarian method (split->hdist after first time) */

    split_remove_duplicates (split->disagree, &(split->n_disagree)); // some elements are equal; this function also qsorts
    split_find_small_disagreement (split);  // could also be one leaf only

    split->spr++;
    if (split->spr > split->cutoff) return split->spr;
    split_remove_small_disagreement (split);

    split_minimize_subtrees (split);
    mismatch = (split->n_g > 0) && (split->n_s > 0); // all edges were in agreement
  }
  return split->spr;
}

void
split_create_agreement_list (splitset split)
{
  int s, g;
  for (g = 0; g < split->n_g; g++) for (s = 0; s < split->n_s; s++)
    if (bipartition_is_equal (split->g_split[g], split->s_split[s])) {
      bipartition_copy (split->agree[split->n_agree++], split->g_split[g]);
      split->n_g--; split_swap_position (split->g_split, g, split->n_g); /* if we don't swap them, we lose ref to "old" value on g_split[] */
      split->n_s--; split_swap_position (split->s_split, s, split->n_s);
      g--; s = split->n_s; /* pretend loop finished, examine again with new values */
    }
  split_remove_agree_edges (split, split->g_split, &(split->n_g));
  split_remove_agree_edges (split, split->s_split, &(split->n_s));
}

void
split_remove_agree_edges (splitset split, bipartition *b, int *nb)
{
  int i, a;
  for (i = 0; i < (*nb); i++) for (a = 0; a < split->n_agree; a++)
    if (bipartition_is_equal (b[i], split->agree[a])) {
      (*nb)--;
      split_swap_position (b, i, (*nb));
      i--;
      a = split->n_agree; /* loop again over new value */
    }
}

void
split_remove_duplicates (bipartition *b, int *nb)
{
  int i, j;
  bipartition pivot;
  if ((*nb) < 2) return; /* only if we have a vector with > 1 element */
  qsort (b, (*nb), sizeof (bipartition), compare_splitset_bipartition_increasing);

  for (i = (*nb) - 1; i >= 1; i--)
    if (bipartition_is_equal (b[i], b[i-1])) {
      pivot = b[i]; /* do not lose a pointer to this element */
      for (j = i; j < (*nb)-1; j++) b[j] = b[j+1];
      b[j] = pivot; /* j = (*nb) - 1, which will become obsolete through next line -->  (*nb)-- */
      (*nb)--;
    }
}

void
split_compress_agreement (splitset split)
{
  int i, j, pair[2];

  for (i = 0; i < split->n_agree; i++) if (split->agree[i]->n_ones == 2) { /* cherry in common, can be represented by just one leaf */
    bipartition_to_int_vector (split->agree[i], pair, 2);
    if (split->leaf_idx) for (j=0; j <= split->size; j++)
      if (split->leaf_idx[j] == pair[1]) split->leaf_idx[j] = pair[0];
    split_remove_redundant_bit (split, pair[1]);
    split_new_size (split,split->agree[0]->n->bits - 1, false); /* false = do not recalculate every bipartition's last elem */
    /* Replacing the removed bit already maintains n_ones and clears the last
     * position in every vector. Recounting all bits after each cherry is redundant. */
    for (j = 0; j < split->n_agree; j++) { /* minimize subtree size and remove single leaves */
      bipartition_flip_to_smaller_set (split->agree[j]); /* agree only */
      if (split->agree[j]->n_ones < 2) split_swap_position (split->agree, j--, --split->n_agree);
    }
    i = -1; /* redo all iterations, with new info (agree[] will be smaller) */
  }
}

void
split_create_disagreement_list (splitset split)
{
  /* Only matched edges can be prune candidates. Keep their bitsets and
   * store the other pairwise comparisons as integer assignment costs. */
  int g, s;
  int n = MAX (split->n_g, split->n_s);
  if (!split->h || split->h->size < n) {
    del_hungarian (split->h);
    split->h = new_hungarian (n);
  }
  /* Real cells are overwritten below. Only the padding needs clearing. */
  for (g=0; g < n; g++) for (s=(g < split->n_g ? split->n_s : 0); s < n; s++)
    hungarian_update_cost (split->h, g, s, 0);
  for (g=0; g < split->n_g; g++) for (s=0; s < split->n_s; s++) {
    bipartition_XOR (split->prune, split->g_split[g], split->s_split[s], true);
    int ones = split->prune->n_ones;
    hungarian_update_cost (split->h, g, s, MIN (ones, split->prune->n->bits - ones));
  }
}

void
split_disagreement_assign_match (splitset split)
{
  int g, max_n = MAX (split->n_g, split->n_s);
  hungarian_solve (split->h, max_n);
  if (split->h->final_cost < 0) biomc2_error ("tree split assignment failed");
  /* now split->h->col_mate will have the pairs */
  /* if we do the matching below it becomes much faster, but we may miss the best prune subtrees in a few cases (do not compromise the algo) */
  split->n_disagree = 0;
  for (g = 0; g < max_n; g++) if ((g < split->n_g) && ( split->h->col_mate[g] < split->n_s)) { /* some matchings might be to dummy edges */
    bipartition_XOR (split->disagree[split->n_disagree], split->g_split[g], split->s_split[split->h->col_mate[g]], true); /* true means to calculate n_ones */
    bipartition_flip_to_smaller_set (split->disagree[split->n_disagree++]);

  }
  if (split->match) { split->hdist_reduced = split->h->final_cost+split->h->initial_cost; split->match = false; }
}

void
split_find_small_disagreement (splitset split)
{
  bipartition dis;
  int a, d;

  bipartition_copy (split->prune, split->disagree[0]); /* smallest, in case we don't find a better one in loop below */
  if (split->prune->n_ones < 2) return;

  dis = new_bipartition_from_bipsize (split->disagree[0]->n);
  for (d = 0; d < split->n_disagree; d++) for (a = 0; a < split->n_agree; a++) {
    if ((split->disagree[d]->n_ones == split->agree[a]->n_ones) ||
        (split->disagree[d]->n_ones == (split->agree[a]->n->bits - split->agree[a]->n_ones))) {
      bipartition_XOR (dis, split->disagree[d], split->agree[a], true);
      if      (!dis->n_ones)               { bipartition_copy (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }
      else if (dis->n_ones == dis->n->bits) { bipartition_NOT (split->prune, split->disagree[d]); d = split->n_disagree; a = split->n_agree; }

    }
  }
  /* check if prune nodes are all on same side of a tree or if they are actually two SPRs (one from each tree) */
  for (d = 0; d < split->n_g; d++) {
    if (!bipartition_contains_bits (split->g_split[d], split->prune)) {
      bipartition_NOT (dis, split->g_split[d]);
      if (!bipartition_contains_bits (dis, split->prune)) { split->spr_extra++; d = split->n_g; }
    }
  }
  del_bipartition (dis);
}

void
split_remove_small_disagreement (splitset split)
{
  int *index, i, j = split->prune->n_ones - 1, k = 0, size = split->agree[0]->n->bits;

  index = (int*) biomc2_malloc (split->prune->n_ones * sizeof (int));
  bipartition_to_int_vector (split->prune, index, split->prune->n_ones);

  if (split->leaf_idx) for (i=0; i <= split->size; i++) {
    if (split->leaf_idx[i] >= 0 && bipartition_is_bit_set (split->prune, split->leaf_idx[i])) {
      split->leaf_recomb[i] = 1;
      split->leaf_idx[i] = -1;
    }
  }
  for (i = size - 1; i >= (size - split->prune->n_ones); i--) {
    if (k > j || index[k] >= (size - split->prune->n_ones)) i = -1;
    else {
      if (i == index[j]) j--;
      else split_replace_bit (split, index[k++], i);
    }
  }

  split_new_size (split,size - split->prune->n_ones, true);
  if (index) free (index);
}

void
split_minimize_subtrees (splitset split)
{
  int i;

  for (i = 0; i < split->n_s; i++) {
    bipartition_flip_to_smaller_set (split->s_split[i]);
    if (split->s_split[i]->n_ones < 2) { split->n_s--; split_swap_position (split->s_split, i, split->n_s); i--; }
  }
  for (i = 0; i < split->n_g; i++) {
    bipartition_flip_to_smaller_set (split->g_split[i]);
    if (split->g_split[i]->n_ones < 2) { split->n_g--; split_swap_position (split->g_split, i, split->n_g); i--; }
  }
  for (i = 0; i < split->n_agree; i++) {
    bipartition_flip_to_smaller_set (split->agree[i]);
    if (split->agree[i]->n_ones < 2) { split->n_agree--; split_swap_position (split->agree, i, split->n_agree); i--; }
  }
}

void
split_remove_redundant_bit (splitset split, int id)
{
  int  last = split->agree[0]->n->bits-1;
  split_replace_bit (split, id, last);
}

void
split_replace_bit (splitset split, int to, int from)
{
  if (from < to) return;
  if (split->leaf_idx) for (int i=0; i <= split->size; i++)
    if (split->leaf_idx[i] == from) split->leaf_idx[i] = to;
  /*not needed for disagree[] */
  bipartition_replace_bit_in_vector (split->agree,   split->n_agree, to, from, true);
  bipartition_replace_bit_in_vector (split->g_split, split->n_g,     to, from, true);
  bipartition_replace_bit_in_vector (split->s_split, split->n_s,     to, from, true);
}

void
split_new_size (splitset split, int size, bool update_bipartitions)
{
  bipsize_resize (split->g_split[0]->n, size);
  bipsize_resize (split->s_split[0]->n, size);
  bipsize_resize (split->agree[0]->n, size);
  bipsize_resize (split->disagree[0]->n, size);
  if (update_bipartitions) {
    bipartition_resize_vector (split->g_split, split->n_g);
    bipartition_resize_vector (split->s_split, split->n_s);
    bipartition_resize_vector (split->agree, split->n_agree);
  }
}

void
split_swap_position (bipartition *b, int i1, int i2)
{
  bipartition pivot = b[i1];
  b[i1] = b[i2];
  b[i2] = pivot;
}
