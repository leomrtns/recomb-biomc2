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


#include "topology.h"

const char *
tree_metric_name (tree_metric metric)
{
  static const char *names[] = {"dspr", "rf", "hdist"};
  return names[metric];
}

tree_metric
parse_tree_metric (const char *name)
{
  for (int i=METRIC_DSPR; i <= METRIC_HDIST; i++)
    if (!strcasecmp (name, tree_metric_name (i))) return i;
  biomc2_error ("unknown tree metric '%s' (expected dspr, rf or hdist)", name);
  return METRIC_DSPR;
}

split_space
new_split_space (int n)
{
  split_space split = biomc2_malloc (sizeof *split);
  split->inodes = n - 3;
  split->n_leaves = n;
  split->metric = METRIC_DSPR;
  split->work = new_splitset_dSPR (n);
  split->leaf_idx = biomc2_malloc (n * sizeof (int));
  split->leaf_recomb = biomc2_malloc (n * sizeof (int));
  split->record_leaves = false;
  return split;
}

void
del_split_space (split_space split)
{
  del_splitset (split->work);
  free (split->leaf_idx);
  free (split->leaf_recomb);
  free (split);
}

void
setup_split_store_recomb_info (split_space split)
{
  split->record_leaves = true;
  for (int i=0; i < split->n_leaves; i++) {
    split->leaf_idx[i] = i;
    split->leaf_recomb[i] = 0;
  }
}

void
setup_split_dont_store_recomb_info (split_space split)
{
  split->record_leaves = false;
}

static void
update_clade (topol_node node)
{
  if (!node->internal) return;
  update_clade (node->left);
  update_clade (node->right);
  bipartition_OR (node->split, node->left->split, node->right->split, false);
}

void
create_topology_bipartition_l (topology tree)
{
  update_clade (tree->root->left);
}

void
order_topology_by_bipartition (topol_node node)
{
  if (!node->internal) return;
  order_topology_by_bipartition (node->left);
  order_topology_by_bipartition (node->right);
  if (bipartition_is_larger (node->left->split, node->right->split)) {
    topol_node tmp = node->left;
    node->left = node->right;
    node->right = tmp;
  }
}

void
euler_semitour_mapping (topol_node node, int *tour, int *n)
{
  /* A fixed child order and internal marker encode shape independently of node IDs. */
  tour[(*n)++] = node->internal ? -1 : node->id;
  if (!node->internal) return;
  euler_semitour_mapping (node->left, tour, n);
  euler_semitour_mapping (node->right, tour, n);
}

static void
prepare_splits (topology tree, bipartition *out)
{
  int n = 0;
  create_topology_bipartition_l (tree);
  bipsize_resize (out[0]->n, tree->nleaves);
  for (int i=tree->nleaves; i < tree->nnodes; i++) {
    topol_node node = tree->nodelist[i];
    if (node == tree->root->left) continue;
    bipartition_copy (out[n], node->split);
    bipartition_flip_to_smaller_set (out[n++]);
  }
  qsort (out, n, sizeof *out, compare_splitset_bipartition_increasing);
}

bool
dSPR_is_zero_l (topology a, topology b, split_space split)
{
  splitset w = split->work;
  prepare_splits (a, w->g_split);
  prepare_splits (b, w->s_split);
  w->n_g = w->n_s = split->inodes;
  for (int i=0; i < w->n_g; i++)
    if (!bipartition_is_equal (w->g_split[i], w->s_split[i])) return false;
  return true;
}

int
dSPR_topology_l (topology a, topology b, split_space split, int *max)
{
  splitset w = split->work;
  int i = 0, j = 0, common = 0;
  split->rfdistance = 0;
  if (dSPR_is_zero_l (a, b, split)) return 0;
  if (!*max) return 1;
  /* The splits are already sorted for equality. Count shared edges in one pass. */
  while (i < w->n_g && j < w->n_s) {
    int cmp = compare_splitset_bipartition_increasing (&w->g_split[i], &w->s_split[j]);
    if (cmp <= 0) i++;
    if (cmp >= 0) j++;
    if (!cmp) common++;
  }
  split->rfdistance = w->n_g - common; /* historical half-RF column in biomc2.distance */
  if (split->metric == METRIC_RF) return 2 * split->rfdistance;
  if (split->metric == METRIC_HDIST) return rf_hdist_topology_lowlevel (w, false);
  if (!split->record_leaves && split->rfdistance == 1) return 1;
  w->leaf_idx = split->record_leaves ? split->leaf_idx : NULL;
  w->leaf_recomb = split->leaf_recomb;
  w->cutoff = *max;
  return dSPR_topology_lowlevel (w);
}

static int
compare_clades (const void *a, const void *b)
{
  topol_node x = *(topol_node const *)a, y = *(topol_node const *)b;
  return compare_splitset_bipartition_increasing (&x->split, &y->split);
}

static void
mark_lower (topol_node node, topology old)
{
  if (!node->internal) return;
  mark_lower (node->left, old);
  mark_lower (node->right, old);
  topol_node prev = old->nodelist[node->map_id];
  int l = node->left->map_id, r = node->right->map_id;
  node->d_done = node->left->d_done && node->right->d_done &&
    ((prev->left->id == l && prev->right->id == r) || (prev->left->id == r && prev->right->id == l));
}

static void
mark_upper (topol_node node, topology old)
{
  if (!node->internal) return;
  if (node->up->internal) {
    topol_node prev = old->nodelist[node->map_id];
    node->u_done = node->up->u_done && node->sister->d_done &&
      prev->up->id == node->up->map_id && prev->sister->id == node->sister->map_id;
  }
  mark_upper (node->left, old);
  mark_upper (node->right, old);
}

void
update_topology_from_topology_split (topology a, topology b, split_space split)
{
  int n = split->inodes + 1, i = 0, j = 0;
  topol_node x[n], y[n];
  create_topology_bipartition_l (a);
  create_topology_bipartition_l (b);
  for (int k=0; k < n; k++) {
    x[k] = a->nodelist[k + a->nleaves];
    y[k] = b->nodelist[k + b->nleaves];
    x[k]->map_id = y[k]->map_id = -1;
  }
  qsort (x, n, sizeof *x, compare_clades);
  qsort (y, n, sizeof *y, compare_clades);
  while (i < n && j < n) {
    int cmp = compare_clades (&x[i], &y[j]);
    if (!cmp) { x[i]->map_id = y[j]->id; y[j]->map_id = x[i]->id; }
    if (cmp <= 0) i++;
    if (cmp >= 0) j++;
  }
  /* Complete the permutation for likelihood buffers belonging to changed edges. */
  j = 0;
  for (i=0; i < n; i++) if (x[i]->map_id < 0) {
    while (y[j]->map_id >= 0) j++;
    x[i]->map_id = y[j]->id;
    y[j]->map_id = x[i]->id;
  }
  a->root->left->u_done = true;
  mark_lower (a->root->left, b);
  mark_upper (a->root->left, b);
}
