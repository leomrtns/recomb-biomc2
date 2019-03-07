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

/*! \file convergence.c
 *  \brief Chain convergence analysis based on Scaled Regeneration Quantile plots.
 */

#include "convergence.h"

int zero = 0;

/* local function prototypes */

/*! \brief Allocate memory for convergence_struct. */
convergence new_convergence (int n_trees, int n_leaves);
/*! \brief Allocate memory for srq_plot_struct. */
srq_plot new_srq_plot (void);
/*! \brief Free memory allocated by convergence_struct. */
void del_convergence (convergence C);
/*! \brief Updates srq_plot_struct if initial state is visited again. */
void update_srq_status (srq_plot SRQ, int *iter);
/*! \brief Calculate relevant statistics (coef. of variation, correlation, slope) of srq_plot_struct. */
void calculate_srq_coefficients (double *cv, double *slope, double *corr_coef, double *freq, int *n_iter, srq_plot SRQ);


convergence
new_convergence (int n_trees, int n_leaves)
{	
	convergence C;
	int i;

	C =
	(convergence) biomc2_malloc (sizeof (struct convergence_struct));
	C->n_t = n_trees;
	C->t = 
	(topology*) biomc2_malloc ((n_trees+1) * sizeof (topology));
	C->position =
	(int*) biomc2_malloc (n_trees * sizeof (int));

	C->SRQ_t =
	(srq_plot*) biomc2_malloc ((n_trees+1) * sizeof (srq_plot));
	C->SRQ_cop =
	(srq_plot*) biomc2_malloc (n_trees * sizeof (srq_plot));

	for (i=0; i < n_trees; i++) {
		C->t[i] = new_topology (n_leaves);
		C->SRQ_t[i] = new_srq_plot ();
    C->SRQ_cop[i] = new_srq_plot ();
		C->position[i] = -1;
	}
	C->t[i] = new_topology (n_leaves);
	C->SRQ_t[i] = new_srq_plot ();

	C->counter = 0;
	return C;
}

srq_plot
new_srq_plot (void)
{
	srq_plot SRQ;
	SRQ = 
	(srq_plot) biomc2_malloc (sizeof (struct srq_plot_struct));
	SRQ->regeneration = SRQ->visit = 0;
	SRQ->Si = SRQ->Sii = SRQ->Sx = SRQ->Sy = SRQ->Sxx = SRQ->Syy = SRQ->Sxy = 0.;
	return SRQ;
}

void
del_convergence (convergence C)
{
	int i;
	if (C) {
		if (C->t) {
			for (i=C->n_t; i >= 0; i--) del_topology (C->t[i]);
			free (C->t);
		}
		if (C->SRQ_t) {
			for (i=C->n_t; i >= 0; i--) if (C->SRQ_t[i]) free (C->SRQ_t[i]); 
			free (C->SRQ_t);
		}
    if (C->SRQ_cop) {
      for (i=C->n_t - 1; i >= 0; i--) if (C->SRQ_cop[i]) free (C->SRQ_cop[i]); 
      free (C->SRQ_cop);
		}
		if (C->position) free (C->position);
	}
	free (C);
}

convergence
prepare_convergence (chain_data chain)
{
	int i=0;
	topology t;
	convergence C = new_convergence (chain->topol->nCOP + 1, chain->ntax);

	for (t=chain->topol->first; (t); t = t->next) {
		copy_topology_from_topology (C->t[i], t);
		C->position[i] = t->last_segment;
		i++;
	}
	copy_topology_from_topology (C->t[i], chain->topol->first);

	return C;
}

void
update_convergence (chain_data chain, convergence C)
{
	topology t;
  int current = 0;

	C->counter++;

  for (t=chain->topol->first; (t) && (current < C->n_t); ) {
    if (t->last_segment >= C->position[current]) {
      /* topologies */
      if (!(dSPR_topology_l (t, C->t[current], chain->split, &zero)) ) {
        update_srq_status (C->SRQ_t[current], &(C->counter));
      }
      /* break-point locations */
			if ( (t->next) && (t->last_segment == C->position[current]) ) {
        update_srq_status (C->SRQ_cop[current], &(C->counter));
      }
      current++;
    }
    else t = t->next;
  }

	/* break-point location if initial stage had nCOP = zero */
	if ( (C->n_t == 1) && (!chain->topol->nCOP)) {
		update_srq_status (C->SRQ_cop[0], &(C->counter));
	}

	/* first topology */
	if (!(dSPR_topology_l (chain->topol->first, C->t[C->n_t], chain->split, &zero)) ) {
		update_srq_status (C->SRQ_t[C->n_t], &(C->counter));
	}
}

void
update_srq_status (srq_plot SRQ, int *iter)
{
	double epslon, dx, dy, di;
  double interval = (double) ( (*iter) - SRQ->regeneration);
    
  SRQ->regeneration = (*iter);
  SRQ->visit++;
  
	if (!SRQ->visit) {
    SRQ->Sx = (double) (SRQ->visit);
    SRQ->Sy = (double) (*iter);
		SRQ->Si = interval; 
    return;
  }

  epslon = (double)(SRQ->visit - 1)/(double)(SRQ->visit);
  dx = (double) (SRQ->visit) - SRQ->Sx;
  dy = (double) (*iter) - SRQ->Sy;
	di = interval - SRQ->Si;

  SRQ->Sxx += (dx * dx * epslon);
  SRQ->Syy += (dy * dy * epslon);
	SRQ->Sxy += (dx * dy * epslon);
  SRQ->Sii += (di * di * epslon);
  SRQ->Sx  += dx/(double)(SRQ->visit);
	SRQ->Sy  += dy/(double)(SRQ->visit);
	SRQ->Si  += di/(double)(SRQ->visit);
}

void
report_convergence (convergence C)
{
	int i;
	double cv, slope, corr_coef, freq;

	if (C->counter < 10) {
		printf ("Sampling too short\n");
		del_convergence (C);
		return;
	}

	printf ("  Scaled Regeneration Quantile plot information [ Y=BX where Y=Ti/TN and X=i/N ]:\n");
	/* first topology */
	calculate_srq_coefficients (&cv, &slope, &corr_coef, &freq, &(C->counter), C->SRQ_t[C->n_t]);
	printf ("freq %-5.4f | coeff.var %-8.6f | slope %-3.2f | corr.coef %6.5f | topology at segment 0\n", 
					freq, cv, slope, corr_coef);
	/* last topologies of each segment */
  for (i=0; i < C->n_t; i++) {
    calculate_srq_coefficients (&cv, &slope, &corr_coef, &freq, &(C->counter), C->SRQ_t[i]);
    printf ("freq %-5.4f | coeff.var %-8.6f | slope %-3.2f | corr.coef %6.5f | topology at segment %-6d\n", 
            freq, cv, slope, corr_coef, C->position[i]);
  }
  for (i=0; i < (C->n_t - 1); i++) {
    calculate_srq_coefficients (&cv, &slope, &corr_coef, &freq, &(C->counter), C->SRQ_cop[i]);
    printf ("freq %-5.4f | coeff.var %-8.6f | slope %-3.2f | corr.coef %6.5f | break-point %-6d\n", 
						freq, cv, slope, corr_coef, C->position[i]);
  }
	/* case where initial state for SRQ had no break-points */
	if (C->n_t == 1) {
		calculate_srq_coefficients (&cv, &slope, &corr_coef, &freq, &(C->counter), C->SRQ_cop[0]);
		printf ("freq %-5.4f | coeff.var %-8.6f | slope %-3.2f | corr.coef %6.5f | no break-points\n", 
						freq, cv, slope, corr_coef);
	}

  del_convergence (C);
}

void 
calculate_srq_coefficients (double *cv, double *slope, double *corr_coef, double *freq, int *n_iter, srq_plot SRQ)
{
	/* scale values to Ti/TN and i/N.
	 * Sxy, Sxx enter the formulae as n*Sxx and n*Sxy so one term visit[i] cancel out. */
	if (SRQ->visit > 5) {
/*
		SRQ->Sx  /= (double) (SRQ->visit);
		SRQ->Sy  /= (double) (SRQ->regeneration);
    SRQ->Sxx /= (double) (SRQ->visit * SRQ->visit); // n Sxx
		SRQ->Syy /= (double) (SRQ->regeneration * SRQ->regeneration);
    SRQ->Sxy /= (double) (SRQ->visit * SRQ->regeneration); // n Sxy
 */
    (*slope) = SRQ->Sxy/SRQ->Sxx;
    (*corr_coef) = SRQ->Sxy/SRQ->Syy;
		(*corr_coef) *= (*slope);
		(*slope) *= ((double)(SRQ->visit)/(double)(SRQ->regeneration));
    if ((*corr_coef) > 0.) (*corr_coef) = sqrt (*corr_coef);
    else (*corr_coef) = 0.;

		(*cv) = (SRQ->Sii/(SRQ->Si * SRQ->Si));
		(*cv) /= (double)(SRQ->visit * SRQ->visit);

		(*freq) = (double)(SRQ->visit)/(double)(*n_iter);
	}
	else (*freq) = (*slope) = (*corr_coef) = (*cv) = 0.;
}
