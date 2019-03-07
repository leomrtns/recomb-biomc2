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

/*! \file summary.c
 * \brief Summarise posterior distribution from biomc2 
 *
 * Reads files "post.tree" and "post.dist" and report point estimates of recombination. This program can be run when
 * posterior samples are still being drawn (before biomc2 finishes). 
 */

#include "summary.h"

void read_summary_topol_from_file (summary sm, char *treefile);
void read_summary_dist_from_file (summary sm, char *distfile);
void remove_unused_memory_from_summary (summary sm);
void remove_duplicate_trees_from_summary (summary sm);
void remove_duplicate_segments_from_summary (summary sm);
void remove_duplicate_samples_from_summary (summary sm);
void create_bp_location_heap (summary sm);
void create_posterior_recomb_freq (summary sm);
void create_empfreq_tree (summary sm);

summary
read_summary_from_file (char *treefile, char *distfile)
{
  summary sm;

  sm = (summary) biomc2_malloc (sizeof (struct summary_struct));

  read_summary_topol_from_file (sm, treefile);
  read_summary_dist_from_file (sm, distfile);
  
  remove_duplicate_segments_from_summary (sm); 
  remove_duplicate_samples_from_summary (sm); 

  /* auxiliary statistics */
  create_bp_location_heap (sm);
  create_posterior_recomb_freq (sm);
  create_empfreq_tree (sm);
  
  return sm;
}

void
read_summary_topol_from_file (summary sm, char *treefile)
{
  nexus_treespace nex;
  int i, *label_order;

  nex = read_nexus_treespace_file (treefile);
  sm->n_index   = sm->n_compact = nex->ntrees;
  sm->t_compact = (topology*) biomc2_malloc (sm->n_index * sizeof (topology));
  sm->t_index   = (int*) biomc2_malloc (sm->n_index * sizeof (int));
  sm->taxlabel  = (char**) biomc2_malloc (sizeof (char*) * nex->nleaves);
  sm->split     = new_split_space (nex->nleaves);
  
  /* order leaves (id numbers) based on t0->taxlabel_hash; mapping is stored in label_order */
  /*  OBS: this functions are useful when you have several tree files and must reorder the leaf IDs */
  label_order = (int*) biomc2_malloc (nex->nleaves * sizeof (int));
  map_hashtable_order (label_order, nex->taxlabel_hash, nex->taxlabel, nex->nleaves);
  order_nexus_treespace_id (nex, label_order);

  for (i=0; i < nex->nleaves; i++) {
    sm->taxlabel[i] = (char*) biomc2_malloc (sizeof (char) * strlen(nex->taxlabel[i]));
    strcpy (sm->taxlabel[i], nex->taxlabel[i]);
  }

  for (i=0; i < sm->n_index; i++) {
    sm->t_compact[i] = new_topology (nex->nleaves);
    copy_topology_from_nexus_tree (sm->t_compact[i], nex->T[i]);
  }
  remove_duplicate_trees_from_summary (sm);
  if (label_order) free (label_order);
  del_nexus_treespace (nex);
}

void
read_summary_dist_from_file (summary sm, char *distfile)
{
  FILE *file;
  int i,j, i1, i2, prev, k, counter = 0;
  char *line_read = NULL, *first_char;
  size_t linelength = 0;

  file = biomc2_fopen (distfile, "r");

  biomc2_getline (&line_read, &linelength, file);
  if (sscanf (line_read, "Nseg %d Nsamples %d", &(sm->n_segments), &(sm->n_samples_original)) != 2)
    biomc2_error ( "malformed post.dist file (first line should be 'Nseg # Nsamples #')\n");

/* sometimes user asks X iter but program outputs (X+1) - rounding error; and program later will remove unused anyhow... */
  sm->n_samples_original++; 

  sm->breakpoint = (int*) biomc2_malloc (sm->n_segments * sizeof (int));
  sm->coldspot   = (int*) biomc2_malloc (sm->n_segments * sizeof (int));

  sm->n_coldspot = 0;

  biomc2_getline (&line_read, &linelength, file);
  first_char = line_read;
  for (i=0; i < sm->n_segments; i++) {
    if (sscanf (first_char, "%d ",&(sm->breakpoint[i])) != 1) 
      biomc2_error ("could not read potential breakpoint (# %d) location (second line of post.dist file)\n", i+1);
    first_char = strstr (first_char, " "); first_char++; 
  }
  sm->n_sites = sm->breakpoint[sm->n_segments-1];

  sm->LnL = (double*) biomc2_malloc (sm->n_samples_original * sizeof (double));
  sm->prior_rate  = (double*) biomc2_malloc (sm->n_samples_original * sizeof (double));
  sm->prior_kappa = (double*) biomc2_malloc (sm->n_samples_original * sizeof (double));
  sm->nSPR = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));
  sm->nCOP = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));

  /* order is different from usual: each element (row) refers to a segment, being a vector of samples (columns) */
  sm->dSPR  = (int**) biomc2_malloc ((sm->n_segments-1) * sizeof (int*));
  sm->topol = (int**) biomc2_malloc (sm->n_segments * sizeof (int*));
  for (i=0; i < sm->n_segments-1; i++) {
    sm->dSPR[i]  = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));
    sm->topol[i] = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));
  }
  sm->topol[i] = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));

  sm->w_sample  = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));
  for (i=0; i < sm->n_samples_original; i++) sm->w_sample[i] = 1;
  sm->sample_index  = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));
  for (i=0; i < sm->n_samples_original; i++) sm->sample_index[i] = i;
  sm->w_segment = (int*) biomc2_malloc (sm->n_segments * sizeof (int));
  for (i=0; i < sm->n_segments; i++) sm->w_segment[i] = 1;

  sm->n_samples = 0;
  while (biomc2_getline (&line_read, &linelength, file) != -1) {
    if (sscanf (line_read, "%d %lf %d %d %lf %lf ", &i1, &(sm->LnL[sm->n_samples]), &(sm->nSPR[sm->n_samples]), 
                &(sm->nCOP[sm->n_samples]), &(sm->prior_rate[sm->n_samples]), &(sm->prior_kappa[sm->n_samples])) != 6) 
      biomc2_error ("could not read initial info for sample on line %d of post.dist file\n", i+3);
    prev = 0;
    first_char = line_read;
    for (j=0; j < sm->nCOP[sm->n_samples]; j++) {
      first_char = strstr (first_char, "|"); first_char++;
      if (sscanf (first_char, " %d %d ",&i1, &i2) != 2) /* i1 refers to the breakpoint and i2 is the SPR distance */ 
        biomc2_error ("could not read breakpoint (# %d) info for sample on line %d of post.dist file\n", j+1, 
                      sm->n_samples+3);
      for (k = prev; k <i1; k++) { 
        sm->topol[k][sm->n_samples] = sm->t_index[counter];
        sm->dSPR[k][sm->n_samples]  = 0;
      }
      sm->topol[k][sm->n_samples] = sm->t_index[counter];
      sm->dSPR[k][sm->n_samples]  = i2;
      prev = i1+1;
      counter++;
    }
    for (k = prev; k <sm->n_segments-1; k++) { 
      sm->topol[k][sm->n_samples] = sm->t_index[counter];
      sm->dSPR[k][sm->n_samples]  = 0;
    }
    sm->topol[k][sm->n_samples] = sm->t_index[counter];
    counter++;
    sm->n_samples++;
  }
  if (!sm->n_samples) biomc2_error ("No posterior samples yet");
  else if (sm->n_samples < sm->n_samples_original) remove_unused_memory_from_summary (sm);

  fclose (file);
  if (line_read) free (line_read);
}

void
del_summary (summary sm)
{
  int i;
  if (!sm) return;

  if (sm->taxlabel) {
    for (i=sm->t_compact[0]->nleaves-1; i>=0; i--)
      if (sm->taxlabel[i]) free (sm->taxlabel[i]);
    free (sm->taxlabel);
  }
  if (sm->t_compact) {
    for (i=sm->n_compact - 1; i >= 0; i--) {
      if (sm->t_compact[i]) del_topology (sm->t_compact[i]);
    }
    free (sm->t_compact);
  }
  if (sm->dSPR) {
    for (i=sm->n_segments-2; i >= 0; i--) if (sm->dSPR[i]) free (sm->dSPR[i]);
    free (sm->dSPR);
  }
  if (sm->topol) {
    for (i=sm->n_segments-1; i >= 0; i--) if (sm->topol[i]) free (sm->topol[i]);
    free (sm->topol);
  }

  if (sm->t_index)    free (sm->t_index);
  if (sm->breakpoint) free (sm->breakpoint);
  if (sm->prior_rate)  free (sm->prior_rate);
  if (sm->prior_kappa) free (sm->prior_kappa);
  if (sm->LnL)  free (sm->LnL);
  if (sm->nSPR) free (sm->nSPR);
  if (sm->nCOP) free (sm->nCOP);
  if (sm->w_sample)  free (sm->w_sample);
  if (sm->w_segment) free (sm->w_segment);
  if (sm->sample_index) free (sm->sample_index);
  del_split_space (sm->split);

  if (sm->bp_loc) {
    for (i=sm->n_samples-1; i >= 0; i--) if (sm->bp_loc[i]) free (sm->bp_loc[i]);
    free (sm->bp_loc);
  }
  if (sm->recomb_freq) free (sm->recomb_freq);
  del_empfreq (sm->nCOP_freq);
  del_empfreq (sm->nSPR_freq);
  if (sm->tree_freq) {
    for (i=sm->n_segments-1; i >= 0; i--) if (sm->tree_freq[i]) del_empfreq (sm->tree_freq[i]);
    free (sm->tree_freq);
  }

  free (sm);
}

void
remove_unused_memory_from_summary (summary sm)
{
  int i;

  sm->n_samples_original = sm->n_samples;

  sm->LnL = (double*) biomc2_realloc ((double*)sm->LnL, sm->n_samples_original * sizeof (double));
  sm->prior_rate  = (double*) biomc2_realloc ((double*) sm->prior_rate, sm->n_samples_original * sizeof (double));
  sm->prior_kappa = (double*) biomc2_realloc ((double*) sm->prior_kappa, sm->n_samples_original * sizeof (double));
  sm->nSPR = (int*) biomc2_realloc ((int*) sm->nSPR, sm->n_samples_original * sizeof (int));
  sm->nCOP = (int*) biomc2_realloc ((int*) sm->nCOP, sm->n_samples_original * sizeof (int));

  for (i=0; i < sm->n_segments-1; i++) {
    sm->dSPR[i]  = (int*) biomc2_realloc ((int*) sm->dSPR[i],  sm->n_samples_original * sizeof (int));
    sm->topol[i] = (int*) biomc2_realloc ((int*) sm->topol[i], sm->n_samples_original * sizeof (int));
  }
  sm->topol[i] = (int*) biomc2_realloc ((int*) sm->topol[i], sm->n_samples_original * sizeof (int));

  sm->w_sample  = (int*) biomc2_realloc ((int*) sm->w_sample, sm->n_samples_original * sizeof (int));
  for (i=0; i < sm->n_samples_original; i++) sm->w_sample[i] = 1;
  sm->sample_index  = (int*) biomc2_realloc ((int*) sm->sample_index, sm->n_samples_original * sizeof (int));
  for (i=0; i < sm->n_samples_original; i++) sm->sample_index[i] = i;
}

void
remove_duplicate_trees_from_summary (summary sm) 
{
  int i, j, *idx;
  topology pivot;
//  clock_t time0, time1;
//  char *s;
//  time0 = clock ();
  
  idx = (int*) biomc2_malloc (sm->n_index * sizeof (int));

  for (i=0; i < sm->n_index; i++) { sm->t_compact[i]->tfreq = 1; idx[i] = sm->t_index[i] = i; }
  for (i=0; i < sm->n_compact - 1; i++) {
    // time1 = clock ();
    // printf ("%6d\t %.12f\n", i, (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    // time0 = time1;
    for (j=i+1; j < sm->n_compact; j++) {
// DEBUG
//      s = topology_to_string (sm->t_compact[i]);
//      printf ("\n%s\t ", s);
//      s = topology_to_string (sm->t_compact[j]);
//      printf ("%s\n", s);
//      free (s);

      if ((dSPR_is_zero_l (sm->t_compact[i], sm->t_compact[j], sm->split))  ) {
        sm->n_compact--;
        // swap duplicate (j) with last topology (n_compact)
        pivot = sm->t_compact[sm->n_compact];
        sm->t_compact[sm->n_compact] = sm->t_compact[j]; 
        sm->t_compact[j] = pivot;
        sm->t_compact[i]->tfreq++;
        // update index
//DEBUG        printf ("DEBUG: [%8d %8d]\n", idx[i], idx[j]);
        sm->t_index[ idx[j] ] = i;
        idx[j] = idx[sm->n_compact];
        j--;
      }
    }
    /* to avoid case where next idx[i] = prev idx[j] = idx[n_compact] */
    if (sm->t_index[ idx[i] ] > i) sm->t_index[ idx[i] ] = i; 
  }
  /* to avoid case where next idx[i] = prev idx[j] = idx[n_compact] */
  if (sm->t_index[ idx[i] ] > i) sm->t_index[ idx[i] ] = i;

  for (i = sm->n_index - 1; i >= sm->n_compact; i--) if (sm->t_compact[i]) del_topology (sm->t_compact[i]);

  free (idx);
}

void
remove_duplicate_segments_from_summary (summary sm)
{
  int i, j, *pivot, n_segments_original = sm->n_segments;
  int freq[n_segments_original-1];

  for (i=0; i < sm->n_segments-1; i++) freq[i] = 0;

  for (i=0; i < sm->n_segments-1; i++)
    for (j=0; j < sm->n_samples; j++)
      freq[i] += sm->dSPR[i][j];

  for (i = sm->n_segments-2; i >= 0; i--) if (!freq[i]) {
    /* shift segments to the "left" (considering the transposed matrix), removing segment [i+1] */

    /* dSPR[i] is redundant */
    pivot = sm->dSPR[i];
    for (j=i; j < sm->n_segments - 2; j++) sm->dSPR[j] = sm->dSPR[j+1];
    sm->dSPR[j] = pivot; /* swap is necessary to avoid memory leak */

    /* topol[i+1] is redundant (since dSPR[i] = distance between topol[i] and topol[i+1] ) */
    pivot = sm->topol[i+1];
    for (j=i+1; j < sm->n_segments - 1; j++) sm->topol[j] = sm->topol[j+1];
    sm->topol[j] = pivot; /* swap is necessary to avoid memory leak */

    /* breakpoint[i] info goes to coldspot (and breakpoint has info about genome size in last segment) */
    sm->coldspot[sm->n_coldspot++] = sm->breakpoint[i];
    for (j=i; j < sm->n_segments - 1; j++) sm->breakpoint[j] = sm->breakpoint[j+1];

    /* w_segment[i] = w_segment[i+1] + 1 (segment [i+1] is the redundant, and w_segment[i]=1 initially ) */
    sm->w_segment[i] += sm->w_segment[i+1];
    for (j=i+1; j < sm->n_segments - 1; j++) sm->w_segment[j] = sm->w_segment[j+1];

    sm->n_segments--;
  }

  /* If the number of segments is reduced, free the unused memory */
  if (sm->n_segments < n_segments_original) {
    sm->breakpoint = (int*) biomc2_realloc ((int*) sm->breakpoint, sm->n_segments * sizeof (int));
    sm->w_segment  = (int*) biomc2_realloc ((int*) sm->w_segment,  sm->n_segments * sizeof (int));
    sm->coldspot   = (int*) biomc2_realloc ((int*) sm->coldspot,   sm->n_coldspot * sizeof (int));
    for (i = n_segments_original - 2; i > sm->n_segments - 2; i--) if (sm->dSPR[i]) free (sm->dSPR[i]);
    sm->dSPR = (int**) biomc2_realloc ((int**) sm->dSPR, (sm->n_segments-1) * sizeof (int*));
    for (i = n_segments_original - 1; i > sm->n_segments - 1; i--) if (sm->topol[i]) free (sm->topol[i]);
    sm->topol = (int**) biomc2_realloc ((int**) sm->topol, sm->n_segments * sizeof (int*));

    /* Sort coldspot values (from first to last segment) */
    qsort (sm->coldspot, sm->n_coldspot, sizeof (int), compare_int);
  }

}

void
remove_duplicate_samples_from_summary (summary sm)
{ 
  int i, j, k, *idx;
  bool equal;

  idx = (int*) biomc2_malloc (sm->n_samples_original * sizeof (int));

  for (i=0; i < sm->n_samples_original; i++) idx[i] = i;

  for (i=0; i < sm->n_samples - 1; i++) {
    for (j=i+1; j < sm->n_samples; j++) {
      for (equal=true, k=0; (equal == true) && (k < sm->n_segments); k++)
        if (sm->topol[k][i] != sm->topol[k][j]) equal = false;

      if (equal == true) { /* samples i and j are identical */
        sm->n_samples--;
        sm->w_sample[i]++;

        if (j < sm->n_samples) {
          sm->nSPR[j] = sm->nSPR[sm->n_samples];
          sm->nCOP[j] = sm->nCOP[sm->n_samples];
          sm->w_sample[j] = sm->w_sample[sm->n_samples];
          for (k=0; k < sm->n_segments-1; k++) {
            sm->topol[k][j] = sm->topol[k][sm->n_samples];
            sm->dSPR[k][j]  = sm->dSPR[k][sm->n_samples];
          }
          sm->topol[k][j] = sm->topol[k][sm->n_samples];
        }
        // update index vector
        sm->sample_index[ idx[j] ] = i;
        idx[j] = idx[sm->n_samples];

        j--; /* since sample [j] now refers to last sample */
      } // if (equal == true)
    } // for (j) 
    sm->sample_index[ idx[i] ] = i;
  } //for (i)
  sm->sample_index[ idx[i] ] = i;

  /* If the number of samples is reduced, free the unused memory */
  if (sm->n_samples < sm->n_samples_original) {
    sm->nSPR = (int*) biomc2_realloc ((int*) sm->nSPR, sm->n_samples * sizeof (int));
    sm->nCOP = (int*) biomc2_realloc ((int*) sm->nCOP, sm->n_samples * sizeof (int));
    sm->w_sample = (int*) biomc2_realloc ((int*) sm->w_sample, sm->n_samples * sizeof (int));
    for (k=0; k < sm->n_segments-1; k++) {
      sm->topol[k] = (int*) biomc2_realloc ((int*) sm->topol[k], sm->n_samples * sizeof (int));
      sm->dSPR[k]  = (int*) biomc2_realloc ((int*) sm->dSPR[k],  sm->n_samples * sizeof (int));
    }
    sm->topol[k] = (int*) biomc2_realloc ((int*) sm->topol[k], sm->n_samples * sizeof (int));
  }
  free (idx);
}

void
create_bp_location_heap (summary sm)
{
  int i,j,k;
  
  sm->bp_loc = (int**) biomc2_malloc (sm->n_samples * sizeof (int*));

  for (i=0; i < sm->n_samples; i++) {
    if (sm->nCOP[i]) {
      sm->bp_loc[i] = (int*) biomc2_malloc ((sm->nCOP[i]+2) * sizeof (int));
      for (k=0, j=0; j < sm->n_segments - 1; j++)
        if(sm->dSPR[j][i]) sm->bp_loc[i][1+k++] = sm->breakpoint[j];
      /* first and last elements are the borders of the genome */
      sm->bp_loc[i][0] = 0; sm->bp_loc[i][sm->nCOP[i]+1] = sm->n_sites;
    }
    else sm->bp_loc[i] = NULL;
  }
}

void
create_posterior_recomb_freq (summary sm)
{
  int i,j;
  sm->nCOP_freq = new_empfreq_from_int_weighted (sm->nCOP, sm->n_samples, sm->w_sample);
  sm->nSPR_freq = new_empfreq_from_int_weighted (sm->nSPR, sm->n_samples, sm->w_sample);

  sm->recomb_freq = (int*) biomc2_malloc ((sm->n_segments - 1) * sizeof (int));
  for (i=0; i < sm->n_segments-1; i++) sm->recomb_freq[i] = 0;
  sm->recomb_freq_sum = sm->recomb_freq_max = 0;

  for (j=0; j < sm->n_samples; j++) 
    for (i=0; i < sm->n_segments-1; i++) if (sm->dSPR[i][j]) {
      sm->recomb_freq[i] += sm->w_sample[j];
      sm->recomb_freq_sum += sm->w_sample[j];
      if (sm->recomb_freq[i] > sm->recomb_freq_max) 
        sm->recomb_freq_max = sm->recomb_freq[i];
    }
}

void
create_empfreq_tree (summary sm)
{
  int i;

  sm->tree_freq = (empfreq*) biomc2_malloc (sm->n_segments * sizeof (empfreq));
  for (i=0; i < sm->n_segments; i++) {
    sm->tree_freq[i] = new_empfreq_from_int_weighted (sm->topol[i], sm->n_samples, sm->w_sample);
// DEBUG
//    printf ("%d [%d] %d\t %d\t .  %d\n",sm->tree_freq[i]->i[0].idx,sm->breakpoint[i], sm->tree_freq[i]->i[0].freq,
//            sm->tree_freq[i]->i[1].idx, sm->tree_freq[i]->i[1].freq); 
  }
//  printf("\n");
}

mosaic
new_mosaic (int n_breakpoints)
{
  mosaic m;
  int i;

  m = (mosaic) biomc2_malloc (sizeof (struct mosaic_struct));
  m->n = n_breakpoints;

  m->loc_idx  = (int*) biomc2_malloc ((n_breakpoints) * sizeof (int));
  m->cd1      = (int*) biomc2_malloc ((n_breakpoints) * sizeof (int));
  m->cd2      = (int*) biomc2_malloc ((n_breakpoints) * sizeof (int));
  m->tree = (empfreq*) biomc2_malloc ((n_breakpoints + 1) * sizeof (empfreq));
  for (i=0; i<=n_breakpoints; i++) m->tree[i] = NULL;

  return m;
}

void
del_mosaic (mosaic m)
{
  int i;

  if (m) {
    if (m->loc_idx) free (m->loc_idx);
    if (m->cd1)     free (m->cd1);
    if (m->cd2)     free (m->cd2);
    if (m->tree) {
      for (i=m->n; i >= 0; i--) del_empfreq (m->tree[i]);
      free (m->tree);
    }
    free (m);
  }
}
