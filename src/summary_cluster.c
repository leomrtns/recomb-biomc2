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

/*! \file summary_cluster.c
 * \brief Calculates distance between breakpoint structures, and finds the centroid structure
 *
 * In the future these functions may find clusters of mosaics and their centroids (which means, using hierarchical
 * clustering algo).
 * 
 */

#include "summary.h"

int pair_distance_bp_loc (int i1, int i2, summary sm);
int single_distance_bp_loc (int *v1, int n1, int *v2, int n2);
void create_mosaic_topology_distribution (mosaic m, summary sm);

dist_matrix
new_distance_matrix (int n_elements, int *weight)
{
  int i;
  dist_matrix m;

  m = (dist_matrix) biomc2_malloc (sizeof (struct dist_matrix_struct));
  m->n = n_elements;
  /* just a pointer, no memory allocation */
  m->w = weight;
  m->d = (int**) biomc2_malloc (n_elements * sizeof (int*));
  for (i=0; i < n_elements; i++) { 
    m->d[i] = (int*) biomc2_malloc (n_elements * sizeof (int));
    m->d[i][i] = 0; /* by definition of distance ;) */
  }

  return m;
}

void
del_distance_matrix (dist_matrix m)
{
  int i;

  if (!m) return;

  if (m->d) {
    for (i= m->n-1; i >= 0; i--) if (m->d[i]) free (m->d[i]);
    free (m->d);
  }
  free (m);
}

void
calculate_distance_matrix_bp_loc (summary sm)
{
  /* At the moment this matrix is not used (see new_mosaic_from_centroid() ) */
  int i, j, tmp = 0, total_dist[sm->n_samples];
  dist_matrix m;
//  empfreq e;

  m = new_distance_matrix (sm->n_samples, sm->w_sample);

  for (i=0; i < sm->n_samples; i++) total_dist[i] = 0;

  for (i=0; i < sm->n_samples-1; i++) {
    for (j=i+1; j < sm->n_samples; j++) {
      m->d[i][j] = m->d[j][i] = pair_distance_bp_loc (i, j, sm);
      total_dist[i] += (m->w[j] * m->d[i][j]);
      total_dist[j] += (m->w[i] * m->d[i][j]);
    }
    if (tmp < total_dist[i]) tmp = total_dist[i]; /* we need a big number */
  }

  /* DEBUG (testing) */
  for (i=0; i < sm->n_samples; i++) if (total_dist[i] < tmp) tmp = total_dist[i];
  for (i=0; i < sm->n_samples; i++) if (total_dist[i] == tmp) {
    for (j=1; j <= sm->nCOP[i]; j++) printf ("%5d ", sm->bp_loc[i][j]);
    printf ("  | total distance: %d\n", tmp);
  }

//  e = new_empfreq_from_int_weighted (total_dist, m->n, m->w);
//  plot_empirical_frequency (e, "empfreq_mosaic_cluster", "mosaic distance distribution", 400., 200., 15.);

//  del_empfreq (e);
  del_distance_matrix (m);
}

int 
pair_distance_bp_loc (int i1, int i2, summary sm)
{
  int i, distance = 0, d1, d2;

  if ((sm->bp_loc[i1] == NULL) && (sm->bp_loc[i2] == NULL)) return 0;
  
  if (sm->bp_loc[i1] == NULL) {
    for (i=1; i < sm->nCOP[i2]+1; i++) {
      d1 = sm->bp_loc[i2][i]; d2 = sm->n_sites - sm->bp_loc[i2][i];
      distance += MIN(d1,d2);
    }
    return distance;
  }
  
  if (sm->bp_loc[i2] == NULL) {
    for (i=1; i < sm->nCOP[i1]+1; i++) {
      d1 = sm->bp_loc[i1][i]; d2 = sm->n_sites - sm->bp_loc[i1][i];
      distance += MIN(d1,d2);
    }
    return distance;
  }

  distance  = single_distance_bp_loc (sm->bp_loc[i1], sm->nCOP[i1], sm->bp_loc[i2], sm->nCOP[i2]);
  distance += single_distance_bp_loc (sm->bp_loc[i2], sm->nCOP[i2], sm->bp_loc[i1], sm->nCOP[i1]);

  return distance;
}

int
single_distance_bp_loc (int *v1, int n1, int *v2, int n2)
{
  int i, j = 1, d1, d2, distance = 0;

  /* v2[0] = 0 and v2[n2+1] = sm->n_sites, so we have two "artificial" breakpoints */
  for (i=1; i < n1+1; i++) {
    d1 = v1[i] - v2[j-1]; d2 = MOD(d2);
    d2 = v1[i] - v2[j]; d2 = MOD(d2);
    while ((d2 < d1) && (j <= n2)) { 
      d1 = d2; 
      d2 = v1[i] - v2[++j]; d2 = MOD(d2);
    }
    if (j > n2) { // reached the border (v2[j] == sm->n_sites)
      distance += MIN (d1, d2);
    }
    else { // d2 became larger than d1  
      distance += d1; 
    }
  }

  return distance;
}

mosaic
new_mosaic_from_centroid (summary sm)
{
  int i, j, dist, dist_sum = 0, min_dist = 0;
  int idx_best, nSPR_best;
  mosaic m;
  empfreq e;
  
  e = new_empfreq (sm->n_samples);

  for (i=0; i < sm->n_samples-1; i++) {
    for (j=i+1; j < sm->n_samples; j++) {
      dist = pair_distance_bp_loc (i, j, sm);
      e->i[i].freq += (sm->w_sample[j] * dist);
      e->i[j].freq += (sm->w_sample[i] * dist);
    } 
    /* the particular value of min_dist don't have any importance, we just need a big number. OTOH, the empfreq->freq
     * values can become really large (can we have overflow?). */
    if (min_dist < e->i[i].freq) min_dist = e->i[i].freq;
  }

  nSPR_best = min_dist = 100 * min_dist; /* just a big number as well */
  for (i=0; i < sm->n_samples; i++) {
    dist_sum += e->i[i].freq;
    if (e->i[i].freq < min_dist) min_dist = e->i[i].freq;
  }

  /* If there is more than one (distinct) ensemble centroid, we pick up the one with smallest overall SPR
   * distance (not only breakpoints). This should work since samples belonging to the centroid usually have identical
   * mosaic structures. Another solution would be to average over all minimum-distance samples. The _real_ solution 
   * would be to find the centroid by (least squares) minimization (since the _real_ centroid may not represented by 
   * any sample. */
  for (i=0; i < sm->n_samples; i++) 
    if ((e->i[i].freq == min_dist) && (sm->nSPR[i] < nSPR_best)) 
     { nSPR_best = sm->nSPR[i]; idx_best = i; }

  m = new_mosaic (sm->nCOP[idx_best]); j = 0;
  /* bp_loc has the site position, but here we need the index so that we can map to the topology vector */
  for (i=0; i < sm->n_segments-1; i++) if (sm->dSPR[i][idx_best]) m->loc_idx[j++] = i;

  /* credible interval based on the 95% samples closest to the centroid */
//  sort_empfreq_increasing (e);

  create_mosaic_topology_distribution (m, sm);

  del_empfreq (e);
  return m;
}

mosaic
new_mosaic_from_quantile (summary sm, double credible_interval)
{
  int i, i1=0, i2=0, i3=0;
  double partial_sum, step_sum, q1_sum, q2_sum, q3_sum; 
  mosaic m = new_mosaic (sm->nCOP_freq->i[0].idx);

  step_sum = (double)(sm->recomb_freq_sum)/(double)(sm->nCOP_freq->i[0].idx);
  q1_sum = ((1.-credible_interval)/2.) * step_sum;
  q2_sum = 0.5   * step_sum; 
  q3_sum = (1. - (1.-credible_interval)/2.) * step_sum;
  /* calculate centered credible intervals and median per breakpoint region */
  partial_sum = sm->recomb_freq[0];
  for(i=0; i < sm->n_segments-2; i++) {
    if (partial_sum > q1_sum) { m->cd1[i1++]     = i; q1_sum += step_sum; }
    if (partial_sum > q2_sum) { m->loc_idx[i2++] = i; q2_sum += step_sum; }
    if (partial_sum > q3_sum) { m->cd2[i3++]     = i; q3_sum += step_sum; }
    partial_sum += sm->recomb_freq[i+1];
  }

  /* just in case (if we have some weird roundoff error or I made some mistake) */
  if (i1 < sm->nCOP_freq->i[0].idx) for (i=i1; i<sm->nCOP_freq->i[0].idx; i++) m->cd1[i] = sm->recomb_freq[sm->n_segments-2];
  if (i2 < sm->nCOP_freq->i[0].idx) for (i=i2; i<sm->nCOP_freq->i[0].idx; i++) m->loc_idx[i] = sm->recomb_freq[sm->n_segments-2];
  if (i3 < sm->nCOP_freq->i[0].idx) for (i=i3; i<sm->nCOP_freq->i[0].idx; i++) m->cd2[i] = sm->recomb_freq[sm->n_segments-2];
  
  create_mosaic_topology_distribution (m, sm);
  return m;
}


void
create_mosaic_topology_distribution (mosaic m, summary sm)
{
  int i, j, k, tfreq[sm->n_compact], first = 0, last, idx[sm->n_compact]; 

  for (i=0; i < sm->n_compact; i++) idx[i] = i;

  for (i=0; i <= m->n; i++) {
    for (k=0; k < sm->n_compact; k++) tfreq[k] = 0;

    if (i < m->n) last = m->loc_idx[i];
    else last = sm->n_segments-1;

    for (j=first; j <= last; j++) {
      for (k=0; k < sm->tree_freq[j]->n; k++) tfreq[ sm->tree_freq[j]->i[k].idx ] += sm->tree_freq[j]->i[k].freq;
      // tfreq[ sm->tree_freq[j]->i[0].idx ] += sm->tree_freq[j]->i[0].freq;
      // if (j==first) printf (".");
      // printf ("%d\t%d\t%d\n", j, sm->tree_freq[j]->i[0].idx, sm->tree_freq[j]->i[0].freq);
    }

    m->tree[i] = new_empfreq_from_int_weighted (idx, sm->n_compact, tfreq);
    first = last + 1;
  }
}
