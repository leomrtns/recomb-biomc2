/*! 
 * Copyright (C) 2006	Leonardo de Oliveira Martins
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

/*! \file summary.h 
 *  \brief Header file for summary.c  
 */

#ifndef _biomc2_summary_h
#define _biomc2_summary_h

#include <cairo.h>
#include <cairo-pdf.h>

#include "empirical_frequency.h"

typedef struct summary_struct* summary;
typedef struct dist_matrix_struct* dist_matrix;
typedef struct mosaic_struct* mosaic;

struct summary_struct 
{
	/*  From post.dist file */

  /*! \brief Number of segments. */
  int n_segments;
  /*! \brief Number of samples (after removal of duplicate samples). */
  int n_samples;
	/*! \brief Genome size (position of last site, reported in post.dist as last breakpoint location). */
	int n_sites;
  /*! \brief Number of original samples (for continuous variables). */
  int n_samples_original;
  /*! \brief Location of possible breakpoints (position of last site of segment). */
  int *breakpoint;
	/*! \brief Location of cold spots (breakpoints with \f$\hat{d}_{SPR}=0\f$ for all samples. */
 	int *coldspot;
	/*! \brief Number of cold spots. */
	int n_coldspot;
  /*! \brief Sample likelihood. */
  double *LnL;
  /*! \brief Sample total SPR distance. */
  int *nSPR;
  /*! \brief Sample total number of breakpoints. */
  int *nCOP;
  /*! \brief Sample prior \f$\mu_0\f$ for substitution rate. */
  double *prior_rate;
  /*! \brief Sample prior \f$\kappa_0\f$ for ti/tv ratio. */
  double *prior_kappa;
  /*! \brief Matrix of SPR distances per sample per segment. */
  int **dSPR;
  /*! \brief Matrix of topology indexes per sample per segment. */
  int **topol;
	/*! \brief Weight of each segment (after identical neighboring segments are removed). */
	int *w_segment;
	/*! \brief Weight of each sample, after removing duplicates. */
	int *w_sample;
  /*! \brief Index of original samples (can be used for convergence purposes, where order is important). */
  int *sample_index;

	/* from post.tree file */

  /*! \brief Number of indexed topologies (total sampled topologies). */
  int n_index;
	/*! \brief Index relating topologies in original sample to summary_struct::t_compact. */
	int *t_index;
  /*! \brief Number of distinct (compacted by frequency) topologies. */
  int n_compact;
	/*! \brief Vector of distinct topologies. */
	topology *t_compact;
	/*! \brief split bipartition information (to calculate topology distance). */
	split_space split;
  /*! \brief labels for tree leaves. */
  char **taxlabel;
  
  /* Descriptive statistics and auxiliary vectors */

  /*! \brief Matrix of breakpoint locations per sample (variable size). */
  int **bp_loc;
  /*! \brief Distribution of breakpoints. */
  empfreq nCOP_freq;
  /*! \brief Distribution of SPR distances. */
  empfreq nSPR_freq;
  /*! \brief posterior frequency (counts) of breakpoints per segment. */
  int *recomb_freq;
  /*! \brief total sum of posterior frequencies (equal to the total number of observed breakpoints). */
  int recomb_freq_sum;
  /*! \brief Count of the most frequent breakpoint. */
  int recomb_freq_max;
  /*! \brief Vector of tree distribution per segment. */
  empfreq *tree_freq;
};

struct dist_matrix_struct
{
  /*! \brief Distance matrix (between mosaics or topology vectors). */
  int **d;
  /*! \brief Number of elements in matrix. */
  int n;
  /*! \brief Weight vector (multiplicity of element). */
  int *w;
};

struct mosaic_struct
{
  /*! \brief Location (index) of breakpoints of the mosaic structure. */
  int *loc_idx;
  /*! \brief Credibility interval (2.5% percentile) of the mosaic structure. */
  int *cd1;
  /*! \brief Credibility interval (97.5% percentile) of the mosaic structure. */
  int *cd2;
  /*! \brief Topology distribution per region. */
  empfreq *tree;
  /*! \brief Size of the mosaic structure (number of breakpoints) given by the quantiles. */
  int n;
};

/* Global function prototypes */

/* from summary.c */
summary read_summary_from_file (char *treefile, char *distfile);
void del_summary (summary sm);

mosaic new_mosaic (int n_breakpoints);
void del_mosaic (mosaic m);

/* from summary_plot.c */
void plot_post_recomb_freq (summary sm);
void plot_empirical_frequency (empfreq e, char *filename, char *title, double width, double height, double border);


/* from summary_cluster.c */
dist_matrix new_distance_matrix (int n_elements, int *weight);
void del_distance_matrix (dist_matrix m);

void calculate_distance_matrix_bp_loc (summary sm);
mosaic new_mosaic_from_centroid (summary sm);
mosaic new_mosaic_from_quantile (summary sm, double credible_interval);

#endif
