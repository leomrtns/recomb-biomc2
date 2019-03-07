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

/*! \file update_chain.h 
 *
 *  \brief Header file for update_chain_params.c and update_chain_topology.c
 */

#ifndef _biomc2_update_chain_h
#define _biomc2_update_chain_h

#include "chain_data.h"

#define WEIGHT 1.

/* Global function prototypes */

/* update_chain_params.c */

void full_likelihood (double *r, chain_data chain);
void lnTruncate (double *r, double *w, double *l, chain_data chain);
void calculate_moments (double *r1, double *r2, int i, chain_data chain);

void propose_swap_chains (chain_data chain1, chain_data chain2, acceptance_frequency *cold_freq, 
													acceptance_frequency *heat_freq);

void proposal_update_lambda  (chain_data chain, int n, acceptance_frequency *afreq);
void proposal_update_penalty (chain_data chain, int n, acceptance_frequency *afreq);
void proposal_update_substitution_rate  (chain_data chain, int n, topology t, acceptance_frequency *afreq);
void proposal_update_substitution_kappa (chain_data chain, int n, topology t, acceptance_frequency *afreq);
void proposal_update_prior_rate  (chain_data chain, acceptance_frequency *afreq);
void proposal_update_prior_kappa (chain_data chain, acceptance_frequency *afreq);

/* update_chain_topology.c */

void mcmc_topol_update_outside (chain_data chain, topology t, acceptance_frequency *afreq);
void proposal_update_topol_birth_forward  (topology t, chain_data chain, acceptance_frequency *afreq, double *log2);
void proposal_update_topol_birth_backward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2);
void proposal_update_topol_death_forward  (topology t, chain_data chain, acceptance_frequency *afreq, double *log2);
void proposal_update_topol_death_backward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2);
void proposal_update_topol_shift_forward  (topology t, chain_data chain, acceptance_frequency *afreq);
void proposal_update_topol_shift_backward (topology t, chain_data chain, acceptance_frequency *afreq);

bool debug_topol (topology t, chain_data chain);

#endif
