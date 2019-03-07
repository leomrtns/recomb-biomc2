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

/*! \file convergence.h 
 *
 *  \brief header file for convergence.c
 */

#ifndef _biomc2_convergence_h
#define _biomc2_convergence_h

#include "update_chain.h"

typedef struct convergence_struct* convergence;
typedef struct srq_plot_struct* srq_plot;

/*! \brief Follow-up of initial state (after warm-up) topologies. */
struct convergence_struct 
{
	/*! \brief Vector with topologies we want to track regeneration times. */
	topology *t;
	/*! \brief Number of topologies. */
	int n_t;
	/*! \brief Vector with positions of topologies (last segment before break-points, 
	 * and first and last topologies). */
	int *position;
	/*! \brief Scaled Regeneration Quantile plot information about topologies. */
	srq_plot *SRQ_t;
	/*! \brief Scaled Regeneration Quantile plot information about break-points. */
	srq_plot *SRQ_cop;
	/*! \brief Iteration counter, to make it independent of MCMC definition of iteration. */
	int counter;
};

/*! \brief Regeneretation times of initial state (after warm-up). */
struct srq_plot_struct
{
	/*! \brief Last time state was visited. */
	int regeneration;
	/*! \brief Number of visits. */
	int visit;
	/*! \brief Sum of visit times. */ 
	double Sx;
	/*! \brief Sum of regeneration times. */ 
	double Sy;
	/*! \brief Sum of squared visit times. */ 
	double Sxx;
	/*! \brief Sum of squared regeneration times. */ 
	double Syy;
	/*! \brief Sum of cross-correlations. */ 
	double Sxy;
	/*! \brief Sum of tour interval (for coefficient of variation calculation). */
	double Si;
	/*! \brief Sum of squared tour interval (for coefficient of variation calculation). */
	double Sii;
};

/* Global function prototypes */

/*! \brief Initialize topologies from convergence_struct based on link_topol_struct. */
convergence prepare_convergence (chain_data chain);
/*! \brief Scan link_topol_struct comparing against convergence_struct. */
void update_convergence (chain_data chain, convergence C);
/*! \brief Print to stdio SRQ information about slope, coef. of variation and correlation coef. */
void report_convergence (convergence C);

#endif
