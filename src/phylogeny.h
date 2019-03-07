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

/*! \file phylogeny.h 
 *
 *  \brief header for phylogeny.c
 */

#ifndef _biomc2_phylogeny_h
#define _biomc2_phylogeny_h

#include "spr.h"

typedef struct phylogeny_struct* phylogeny;
typedef struct node_likelihood_struct* node_likelihood;
typedef struct lk_vector_struct* lk_vector;

/*! \brief Model parameters and likelihood vectors for one segment. */
struct phylogeny_struct
{
	/*! \brief Likelihood vector for each node. */
	node_likelihood *l;
	/*! \brief Vector of pointers to likelihood, used by swap_likelihood().
	 * Necessary not to loose track of old configuration. */ 
	node_likelihood *l_ptr;
  /*! \brief Sequence segment size (in number of distinct patterns, not sites). */
	int nchar;
	/*! \brief Number of taxa. */
	int ntax;
	/*! \brief Number of nodes (internal and leaves). */
	int nnodes;
	/*! \brief Current \f$ ln(L) \f$. Equivalent to phylogeny_struct::likelihood_accepted 
	 * within a mini-sampler. */
	double likelihood_current;
	/*! \brief Proposal \f$ ln(L) \f$. Ultimately subject to acceptance/rejection
	 * by MCMC.*/
	double likelihood_proposal;
	/*! \brief Accepted \f$ ln(L) \f$. */
	double likelihood_accepted;
	/*! \brief Current \f$ \mu_i \f$ (expected substitution rate). */
	double rate;
	/*! \brief Current transition/transversion ratio \f$\kappa_i\f$ for HKY model. */
	double kappa;
	/*! \brief Estimated topology distance \f$d_{SPR}\f$ to next segment. */ 
  int dist;
  /*! \brief Prior \f$\lambda_i\f$ for topology distances (modified truncated Poisson). */
  double lambda;
  /*! \brief Penalty weight \f$w_i\f$ term of modified Poisson (corrects for underdispersion). */ 
  double penalty;
	/*! \brief Frequency of each site pattern (how many sites have the same AGCT pattern). */
	double *site_weight;
  /*! \brief Transition probability matrix. */
	double Q[4][4];
	/*! \brief Equilibrium base distribution shared amongst segments (pointer to chain_data_struct::pi). */
	double *pi;
	double *z1[4], /*!< \brief Left eigenvector for HKY model (pointer to chain_data_struct::z1). */
         *z2[4]; /*!< \brief Right eigenvector for HKY model (pointer to chain_data_struct::z2). */
  /*! \brief Eigenvalues for HKY model (specific to segment since is a function of \f$\mu_i\f$ and \f$\kappa_i\f$). */
	double psi[4];
};

/*! \brief Partial Likelihood information for each node such that no calculation is necessary 
 * if new state is rejected. */
struct node_likelihood_struct
{
  /*! \brief Upstream partial likelihood linked list. */
  lk_vector *u;
  /*! \brief Dounstream partial likelihood linked list. */
  lk_vector *d;
  /*! \brief Likelihood linked list size (1 for leaves and largest minisampler size for internal nodes). */
  int size;
  /*! \brief Upstream partial likelihood of current topology (inside Al-Awadi update, for instance). Proposal topologies are 
   * accessed through node_likelihood_struct::u_current->next. */
	lk_vector u_current, 
  /*! \brief downstream partial likelihood of current topology (inside Al-Awadi update, for instance). Proposal topologies are 
   * accessed through node_likelihood_struct::d_current->next. */
            d_current;
  lk_vector u_accepted, /*!< \brief Upstream partial likelihood of last accepted topology (before update). */
            d_accepted; /*!< \brief Downstream partial likelihood of last accepted topology (before update). */
};

/*! \brief Circular linked list with partial likelihood information for a node. Its size is the largest between 
 * chain_data_struct::n_cycles and chain_data_struct::n_mini. */
struct lk_vector_struct
{
  /*! \brief Partial likelihood values for each pattern and state (A,G,C,T). */
  double **lk;
  /*! \brief Double-linked circular list information (lk_vector_struct::next
   * represents proposal, and accept*/
	lk_vector next, prev;
};

/* Global function prototypes */

phylogeny new_phylogeny (int nchar, int ntax, int ncycles);
void del_phylogeny (phylogeny p);

phylogeny new_phylogeny_from_nexus (nexus_alignment a, double *siteweight, int *idx, int n_idx, int ncycles);
double average_rate_distance_method (phylogeny phy);
void update_equilibrium_frequencies (double *pi, phylogeny phy);

#endif
