/*!
 * Copyright (C) 2006 Leonardo de Oliveira Martins
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

/*! \file run_sampler.c
 *
 *  \brief Run heat and cold chains in serial.
 */

#include "run_sampler.h"

/* local function prototypes */

void print_initial_info(chain_data chain, int interval);
void print_header (void);
void set_frequency_to_zero (acceptance_frequency *afreq);
void print_acceptance_rates (chain_data chain, acceptance_frequency *afreq);
void estimate_time (int seconds, double ntimes);
void mcmc_setup_initial_values (chain_data chain);
void mcmc_set_temperature (chain_data chain1, chain_data chain2, int *i, int *j);
void mcmc_forward_backward (chain_data chain, int *iter, acceptance_frequency *afreq, convergence C);
void mcmc_forward  (chain_data chain, acceptance_frequency *afreq);
void mcmc_backward (chain_data chain, acceptance_frequency *afreq);

void clean_linked_topology (chain_data chain, int *iter);

void sample_to_output_file (chain_data chain_file, chain_data chain, int *j);
void open_output_files (chain_data chain);
void close_output_files (chain_data chain);

void 
run_sampler (char *filename)
{
  int i, j, time_estim = 0, interval;
//  int n_mini_bkp, n_cycles_bkp;
//  double kT_mini_bkp, kT_cycles_bkp;
  clock_t time0, time1;
  acceptance_frequency cold_freq, heat_freq;
  chain_data chain1, chain2;
  convergence C_cold, C_heat;

  chain1 = new_chain_data_from_file (filename);
  chain2 = new_chain_data_from_file (filename);

  setup_acceptance_names ();
  mcmc_setup_initial_values (chain1);
  mcmc_setup_initial_values (chain2);

  time1 = clock ();

  interval = chain1->ngen/chain1->nsamples;
  chain2->nsamples = chain1->nsamples;

  print_initial_info (chain1, interval);

	set_frequency_to_zero ( &cold_freq);
	set_frequency_to_zero ( &heat_freq);

	/* Posterior topologies */
	open_output_files (chain1);

	/* stdout header */
	printf ("WARM-UP STAGE (simulated annealing)\n");

  for (i = 1; i <= 5; i++) {
		/* pre- burn-in to unlink topologies on heated chain */
		for (j = 1; j <= chain1->warmup; j++) { 
			if (!(j%chain1->ratio_spr)) {
				branch_swap = &(spr_apply_random_move_spr);
				chain1->logD = &(chain1->logDspr);
				chain2->logD = &(chain2->logDspr);
			}
			else {
				branch_swap = &(spr_apply_random_move_nni);
				chain1->logD = &(chain1->logDnni);
				chain2->logD = &(chain2->logDnni);
			}
			
			mcmc_set_temperature (chain1, chain2, &i, &j);
			mcmc_forward_backward (chain1, &j, &cold_freq, NULL);
			mcmc_forward_backward (chain2, &j, &heat_freq, NULL);
		}
		
		time0 = time1; time1 = clock (); 

    printf ("\n[warm-up] iteration %5d, %8.4f secs, ", i*(j-1), (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    printf ("cooling scheme (1/temperature): %f - %f \n", chain1->beta_zero, chain1->beta_n);
		print_header ();
    printf ("ch 1: "); print_acceptance_rates (chain1, &cold_freq);
		printf ("ch 2: "); print_acceptance_rates (chain2, &heat_freq);
		set_frequency_to_zero ( &cold_freq);
		set_frequency_to_zero ( &heat_freq);

		time_estim += (int)(time1 - time0);
	} // Warm-up 

/*  
  chain1->n_mini    = chain2->n_mini    = n_mini_bkp;
  chain1->kT_mini   = chain2->kT_mini   = kT_mini_bkp;
  chain1->n_cycles  = chain2->n_cycles  = n_cycles_bkp;
  chain1->kT_cycles = chain2->kT_cycles = kT_cycles_bkp;
	chain1->acceptance_probability = &(acceptance_probability_bayesian);
	chain2->acceptance_probability = &(acceptance_probability_bayesian);
*/

	printf ("\nEstimated time to completion : "); 
  estimate_time ( time_estim, (double)(chain1->ngen + (5 * chain1->burnin))/(double)(5 * chain1->warmup) );
	printf ("\n\n");
  time_estim = 0;
  chain1->kT = chain2->kT = 1.;
	/* Setup SRQ initial states based on independent cold chains */
  printf ("BURN-IN STAGE (no swap between cold and heated chains): ");
  printf ("Sampling initial state for convergence statistics\n");
  for (i = 1; i <= 5; i++) {
    for (j = 1; j <= chain1->burnin; j++) { 
			if (!(j%chain1->ratio_spr)) {
				branch_swap = &(spr_apply_random_move_spr);
				chain1->logD = &(chain1->logDspr);
				chain2->logD = &(chain2->logDspr);
			}
			else {
				branch_swap = &(spr_apply_random_move_nni);
				chain1->logD = &(chain1->logDnni);
				chain2->logD = &(chain2->logDnni);
			}
  
			mcmc_forward_backward (chain1, &j, &cold_freq, NULL);
      mcmc_forward_backward (chain2, &j, &heat_freq, NULL);
    }
    time0 = time1; time1 = clock (); 
    printf ("\n[burn-in] iteration %5d, %8.4f secs\n", i*(j-1), (double)(time1-time0)/(double)CLOCKS_PER_SEC);
    print_header ();
    printf ("ch 1: "); print_acceptance_rates (chain1, &cold_freq);
    printf ("ch 2: "); print_acceptance_rates (chain2, &heat_freq);
    set_frequency_to_zero ( &cold_freq);
    set_frequency_to_zero ( &heat_freq);
    time_estim += (int)(time1 - time0);
  }
  C_cold = prepare_convergence (chain1);
  C_heat = prepare_convergence (chain2);

  printf ("\nEstimated time to completion : "); 
  //	estimate_time ( time_estim, (double)(chain1->ngen)/(double)(10 * chain1->warmup) );
  estimate_time ( time_estim, (double)(chain1->ngen)/(double)(5 * chain1->burnin) );
  printf ("\n\n");

	/* main sampler iteration */
	printf ("\nMAIN SAMPLER (saving prior/posterior parameters/topologies to file)\n");
	chain1->kT = 1.;              /* cold chain */
	chain2->kT = chain2->kT_heat; /* heated chain */

	for (j = 1; j <= chain1->ngen; j++) {
		if (!(j%chain1->ratio_spr)) {
			branch_swap = &(spr_apply_random_move_spr);
			chain1->logD = &(chain1->logDspr);
			chain2->logD = &(chain2->logDspr);
		}
		else {
			branch_swap = &(spr_apply_random_move_nni);
			chain1->logD = &(chain1->logDnni);
			chain2->logD = &(chain2->logDnni);
		}

		if (chain1->kT == 1.) {
			mcmc_forward_backward (chain1, &j, &cold_freq, C_cold);
			mcmc_forward_backward (chain2, &j, &heat_freq, C_heat);
		}
		else {
			mcmc_forward_backward (chain2, &j, &cold_freq, C_cold);
			mcmc_forward_backward (chain1, &j, &heat_freq, C_heat);
		}
		
		if (!(j%chain1->swap_interval)) propose_swap_chains (chain1, chain2, &cold_freq, &heat_freq);
 
    if(!(j%interval)) {
      /* stream IO info is stored only in chain1 */
      if (chain1->kT == 1.) sample_to_output_file (chain1, chain1, &j);
      else sample_to_output_file (chain1, chain2, &j);
    }

    /* show accept_prob 20 times */
		if (!(j%(chain1->ngen/20))) {
			time0 = time1;
			time1 = clock ();
			
			printf ("\n[main sampler] iteration %d, remaining time : ", j);
			estimate_time ( (int)(time1 - time0), (double)((chain1->ngen - j) * 20)/(double)(chain1->ngen) );
			printf ("\n");
			print_header ();

			if (chain1->kT == 1.) {
				printf ("cold: "); print_acceptance_rates (chain1, &cold_freq);
				printf ("heat: "); print_acceptance_rates (chain2, &heat_freq);
			}
			else {
				printf ("cold: "); print_acceptance_rates (chain2, &cold_freq);
				printf ("heat: "); print_acceptance_rates (chain1, &heat_freq);
			}
			set_frequency_to_zero ( &cold_freq);
			set_frequency_to_zero ( &heat_freq);
		}
	}
	printf ("\n [cold chain] "); report_convergence (C_cold);
	printf ("\n [heat chain] "); report_convergence (C_heat);

  close_output_files (chain1);
  del_chain_data (chain1);
  del_chain_data (chain2);
}

void
print_initial_info(chain_data chain, int interval)
{
	printf ("  number of segments = %d\n", chain->n_segments);
  printf ("  ngen = %d burnin = %d(x10) sample_interval (after burnin) = %d warm-up = %d(x10) \n", 
          chain->ngen, chain->burnin, interval, chain->warmup);
	printf ("  period, in number of iterations, between SPR moves (otherwise NNI): %d \n",  chain->ratio_spr);
	printf ("- in Al-Awadhi updates (changing the break-point structure):\n");
	printf ("            temperature = %f x normal, mini-sampler size = %d \n", chain->kT_cycles, chain->n_cycles);
	printf ("- in topology updates (without changing the break-point structure):\n");
	printf ("            temperature = %f x normal, mini-sampler size = %d \n", chain->kT_mini, chain->n_mini);
  printf ("            how many times, per iteration, to update topologies: %d \n", chain->ratio_topol_update);  
	printf ("  heated chain: heat factor 1/kT = %f, swap interval = %d\n\n", chain->kT_heat, chain->swap_interval);
	printf (" the acceptance rates below are (number of accepted)/(number of tries);\n");
	printf (" prop means (number of tries for this move)/(total number of tries for moves)\n\n");
}	

void
print_header (void)
{
	int i;
	printf ("      ");
	for (i=0; i < 3; i++) printf ("%s prop|", a_names[i]);
	for (i=3; i < AFREQSIZE; i++) printf (" %s", a_names[i]);
	printf ("\n");
}

void
set_frequency_to_zero (acceptance_frequency *afreq)
{
	int i;
	for (i=0; i < AFREQSIZE; i++) afreq->propose[i] = afreq->accept[i] = 0;
}	

void
print_acceptance_rates (chain_data chain, acceptance_frequency *afreq)
{
	int i;
  double total = 0., fraction = 0.;

  for (i=0; i<3;i++) total += (double)(afreq->propose[i]);
  if (total < 1.) total=1.;

  for (i=0; i < 3; i++) {	
    if (afreq->propose[i]) fraction = (double)(afreq->accept[i])/(double)(afreq->propose[i]);
    else fraction = 0.;
    printf ("%.5f %.2f|", fraction, (double)(afreq->propose[i])/total);
  } 
  printf (" ");
  for (i=3; i < AFREQSIZE; i++) {	
    if (afreq->propose[i]) fraction = (double)(afreq->accept[i])/(double)(afreq->propose[i]);
    else fraction = 0.;
    printf ("%.5f ", fraction);
  }
  printf (" [%d %d]\n", chain->topol->nSPR, chain->topol->nCOP);
}

void
estimate_time (int seconds, double ntimes)
{
	seconds = (int)(((double)(seconds) * ntimes)/(double)(CLOCKS_PER_SEC));
	if (seconds > 86400) { printf ("%2d d ", seconds/86400); seconds = (seconds%86400); }
	if (seconds > 3600)  { printf ("%2d h ", seconds/3600);  seconds = (seconds%3600); }
	if (seconds > 60)    { printf ("%2d m ", seconds/60);    seconds = (seconds%60); }
	printf ("%2d sec ", seconds); 
}

void
mcmc_setup_initial_values (chain_data chain)
{
  int i;
	double w;

	for (i=0; i < chain->n_segments - 1; i++) {
		w = chain->segment[i]->penalty + WEIGHT;
		lnTruncate (&(chain->lnN[i]), &w, &(chain->segment[i]->lambda), chain);
	}

} 

void
mcmc_set_temperature (chain_data chain1, chain_data chain2, int *i, int *j)
{ 
/*
	chain1->beta_n = chain1->beta_zero + (log ( (double)((((*i)-1)*chain1->warmup) + (*j)) ))/4.;
	chain2->beta_n = exp (chain1->beta_n);
	chain1->beta_n = chain2->beta_n;
 */
  (void) i;
  chain1->beta_n = chain1->beta_zero * log ( (double)( (*j)) + EXP_1 );
  chain1->kT = chain2->kT = chain1->beta_n;
}

void
mcmc_forward_backward (chain_data chain, int *iter, acceptance_frequency *afreq, convergence C)
{
	int i;
	topology t;

	mcmc_forward (chain, afreq);
	clean_linked_topology (chain, iter);

	for (i = 0; i < chain->ratio_topol_update; i++) {
		for (t=chain->topol->first; (t); t = t->next) mcmc_topol_update_outside (chain, t, afreq);
		clean_linked_topology (chain, iter);
	}
	if (C) update_convergence (chain, C);

	mcmc_backward (chain, afreq);
	clean_linked_topology (chain, iter);

	for (i = 0; i < chain->ratio_topol_update; i++) {
		for (t=chain->topol->last; (t); t = t->prev) mcmc_topol_update_outside (chain, t, afreq);
		clean_linked_topology (chain, iter);
	}
	if (C) update_convergence (chain, C);
 

	if (!((*iter)%5)) {
		for (t=chain->topol->first; (t); t = t->next) 
			for (i=t->first_segment; i <= t->last_segment; i++) {
				proposal_update_substitution_rate  (chain, i, t, afreq);
				proposal_update_substitution_kappa (chain, i, t, afreq);
			}
	}

	for (i=0; i < chain->n_segments - 1; i++) {
		proposal_update_lambda  (chain, i, afreq);
		proposal_update_penalty (chain, i, afreq);
	}

	if (!((*iter)%2)) {
		proposal_update_prior_rate  (chain, afreq);
		proposal_update_prior_kappa (chain, afreq);
	}

}

void
mcmc_forward (chain_data chain, acceptance_frequency *afreq)
{
	topology t = chain->topol->first;
	double u, log2 = chain->log2;

	if (t->last_segment > t->first_segment) {
		if (!(t->next)) log2 = 0.;
		proposal_update_topol_birth_forward (t, chain, afreq, &log2);
	}
	t = t->next;

	while (t) {
		log2 = 0.;
		
		if (t->last_segment > t->first_segment) {
			u = biomc2_rng_uniform (random_number);
			if (u < 0.4) {
//				if (!(t->next)) log2 = chain->log2;
				proposal_update_topol_birth_forward (t, chain, afreq, &log2);
			}
			else if (u < 0.8) {
				if (!(t->next)) log2 = -chain->log2;
				proposal_update_topol_death_forward (t, chain, afreq, &log2);
			}
			else {
				proposal_update_topol_shift_forward (t, chain, afreq);
			}
		}
		
		else {
			if (!(t->next)) log2 = -chain->log2;
			proposal_update_topol_death_forward (t, chain, afreq, &log2);
		}

		if ( (t->next) && (!(chain->segment[t->last_segment]->dist)) ) { t = t->next; }
    t = t->next;
	}
}

void
mcmc_backward (chain_data chain, acceptance_frequency *afreq)
{
	topology t = chain->topol->last;
	double u, log2 = chain->log2;

	if (t->last_segment > t->first_segment) {
		if (!(t->prev)) log2 = 0.;
		proposal_update_topol_birth_backward (t, chain, afreq, &log2);
	}
	t = t->prev;

	while (t) {
		log2 = 0.;
		
		if (t->last_segment > t->first_segment) {
			u = biomc2_rng_uniform (random_number);
			if (u < 0.4) {
//				if (!(t->prev)) log2 = chain->log2;
				proposal_update_topol_birth_backward (t, chain, afreq, &log2);
			}
			else if (u < 0.8) {
				if (!(t->prev)) log2 = -chain->log2;
				proposal_update_topol_death_backward (t, chain, afreq, &log2);
			}
			else {
				proposal_update_topol_shift_backward (t, chain, afreq);
			}
		}
	
		else {
			if (!(t->prev)) log2 = -chain->log2;
			proposal_update_topol_death_backward (t, chain, afreq, &log2);
		}

		if ( (t->prev) && (!(chain->segment[t->prev->last_segment]->dist)) ) { t = t->prev; }
		t = t->prev;
	}
}

void
clean_linked_topology (chain_data chain, int *iter)
{
	int r_lt_l = 0;
	topology t, redundant, keep;
  (void) iter;

	chain->topol->nSPR = chain->topol->nCOP = 0;
	for (t = chain->topol->first; (t->next); ) {
		if (chain->segment[t->last_segment]->dist) {
			chain->topol->nCOP++; chain->topol->nSPR += chain->segment[t->last_segment]->dist;
			t = t->next;
		}
		else { /* the two segments will be merged (since dist = 0) */
			/* r_lt_l = right topol larger (more segments) than left topol */
			r_lt_l = t->next->last_segment - t->last_segment + t->first_segment - t->next->first_segment;
			if (r_lt_l > 0) 
			 { redundant = t; keep = t->next; }
			else 
			 { redundant = t->next; keep = t; }
			
			/* swap likelihood nodes */
			copy_topology_from_topology_mapping (redundant, keep, chain->split);
			swap_likelihood (redundant, chain->segment, keep);
/*
			redundant->likelihood_proposal = redundant->likelihood_accepted;
			if (debug_topol (redundant, chain)) { 
				printf ("redundant\n");
			}
			redundant->likelihood_accepted = redundant->likelihood_proposal;
 */
			keep->likelihood_accepted += redundant->likelihood_accepted;
			if (r_lt_l > 0) keep->first_segment = redundant->first_segment;
			else keep->last_segment = redundant->last_segment;
			
			remove_topology (redundant, chain->topol);
			t = keep;
/*		
			keep->likelihood_proposal = keep->likelihood_accepted;
			if (debug_topol (keep, chain)) {
				printf ("keep\n");
			}
			keep->likelihood_accepted = keep->likelihood_proposal;
 */
		}
	}
}

void
sample_to_output_file (chain_data chain_file, chain_data chain, int *j)
{
  int i=0;
  char *s;
  double lnLik = 0.;
	topology t;

	for (t=chain->topol->first; (t); t = t->next) {
		create_topology_bipartition_l (t);
		order_topology_by_bipartition (t->root->left);
		s = topology_to_string (t);
		fprintf (chain_file->treefile, "tree rep%ds%d = %s\n",(*j), i++, s);
		fflush (chain_file->treefile);
		free (s);
	}

  fprintf (chain_file->distfile, "%-8d ", (*j));

  full_likelihood (&lnLik, chain);
	fprintf (chain_file->distfile, "%12f ", lnLik);
	fprintf (chain_file->distfile, "%4d ", chain->topol->nSPR);
	fprintf (chain_file->distfile, "%4d ", chain->topol->nCOP);
  fprintf (chain_file->distfile, "%12f ", chain->prior_rate);
  fprintf (chain_file->distfile, "%12f ", chain->prior_kappa);

	for (t=chain->topol->first; (t->next); t = t->next) { 
		fprintf (chain_file->distfile, "| %4d %2d ", t->last_segment, chain->segment[t->last_segment]->dist);
	}
	fprintf (chain_file->distfile, "\n");
	fflush (chain_file->distfile);

}

void
open_output_files (chain_data chain)
{
	int i;

	/* posterior parameters */
	if (chain->prior) {
		chain->distfile = biomc2_fopen ("prior.dist", "w");
		chain->treefile = biomc2_fopen ("prior.tree", "w");
	} else {
		chain->distfile = biomc2_fopen ("post.dist", "w");
		chain->treefile = biomc2_fopen ("post.tree", "w");
	}

	fprintf (chain->treefile, "#NEXUS\n\nBegin trees;\n Translate\n");
	fprintf (chain->treefile, "\t1  %s", chain->taxlabel[0]);
	for (i=1; i < chain->ntax; i++) 
		fprintf (chain->treefile, ",\n\t%d  %s", i+1, chain->taxlabel[i]);
	fprintf (chain->treefile, "\n;\n");


  /* file Header */ 
  fprintf (chain->distfile, "Nseg %d Nsamples %d \n",chain->n_segments, chain->nsamples);
  for (i=0; i<chain->n_segments-1; i++) fprintf (chain->distfile, "%d ", chain->breakpoint[i]);
	fprintf (chain->distfile, "%d\n", chain->nsites);
}

void
close_output_files (chain_data chain)
{
	fprintf (chain->treefile, "\nEnd;\n");
	fclose (chain->treefile); 
  fclose (chain->distfile);
}
