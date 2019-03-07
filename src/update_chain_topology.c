#include "update_chain.h"

/* local function prototypes */

void poisson_dist_term (double *result, chain_data chain, phylogeny p, int *new_dist);
void poisson_dist_term_awadhi (double *result, chain_data chain, phylogeny p, int *new_dist, double *kT1, double *kT2);
void poisson_dist_term_no_ratio (double *result, chain_data chain, phylogeny p);

int mcmc_topol_update (chain_data chain, topology t);
int mcmc_topol_update_shuffle (chain_data chain, topology t);
void al_awadhi_minisampler (chain_data chain, topology t, acceptance_frequency *afreq);

void
poisson_dist_term (double *result, chain_data chain, phylogeny p, int *new_dist)
{
	(*result) =  ( chain->factorial[p->dist] - chain->factorial[(*new_dist)] );
	(*result)	+= ( log (p->lambda) * (double)((*new_dist) - p->dist) ); 
	(*result)	*= (p->penalty + WEIGHT);
}

void
poisson_dist_term_awadhi (double *result, chain_data chain, phylogeny p, int *new_dist, double *kT1, double *kT2)
{
	(*result)	 = ((log (p->lambda)*(double)(*new_dist)) - chain->factorial[*new_dist]) * ((p->penalty + WEIGHT) * (*kT1));
	(*result)	-= ( ((log (p->lambda)*(double)(p->dist)) - chain->factorial[p->dist]) * ((p->penalty + WEIGHT) * (*kT2)) );
}

void
poisson_dist_term_no_ratio (double *result, chain_data chain, phylogeny p)
{
	(*result)	= ((log (p->lambda) * (double)(p->dist)) - chain->factorial[p->dist]) * (p->penalty + WEIGHT);
}

int
mcmc_topol_update (chain_data chain, topology t)
{
  double sum_l = 0., sum_r = 0., acceptance, transition = 0.;
	int d_lp = -1, d_rp = -1;

  (*branch_swap) (t, false);

	if (t->prev) {
		d_lp = dSPR_topology_l (t, t->prev, chain->split, &(chain->max_dist_split));
		if (d_lp > chain->max_distance) {
			spr_apply_move_at_nodes (t, t->undo_prune, t->undo_regraft, true); return 0;
		}
		else poisson_dist_term (&sum_l, chain, chain->segment[t->prev->last_segment], &d_lp);
	}
	if (t->next) {
		d_rp = dSPR_topology_l (t->next, t, chain->split, &(chain->max_dist_split));
		if (d_rp > chain->max_distance) { 
			spr_apply_move_at_nodes (t, t->undo_prune, t->undo_regraft, true); return 0;
		}
		else poisson_dist_term (&sum_r, chain, chain->segment[t->last_segment], &d_rp);
	}

	chain->ln_likelihood_moved_branches (chain->segment, t);
//	if (debug_topol (t, chain)) printf ("topol_update\n");

	acceptance = chain->kT * ( (t->likelihood_proposal - t->likelihood_current) + sum_l + sum_r);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		accept_likelihood_moved_branches (chain->segment, t);

		if (t->prev) chain->segment[t->prev->last_segment]->dist = d_lp;
		if (t->next) chain->segment[t->last_segment]->dist = d_rp;
		return 1;
	}
	else { 
		spr_apply_move_at_nodes (t, t->undo_prune, t->undo_regraft, true);
		clear_topology (t);
 		return 0;
	}
}

int
mcmc_topol_update_shuffle (chain_data chain, topology t)
{
  double sum_l = 0., sum_r = 0., acceptance, transition = 0.;
	int d_lp = -1, d_rp = -1, i, n_shuffle = biomc2_rng_uniform_int (random_number, 3);

	copy_topology_from_topology (chain->topol->bkp[BKPshuffle_random], t);
	for (i=0; i <= n_shuffle; i++) 
		spr_apply_random_move_spr (chain->topol->bkp[BKPshuffle_random], true);

	if (t->prev) {
		d_lp = dSPR_topology_l (chain->topol->bkp[BKPshuffle_random], t->prev, chain->split, &(chain->max_dist_split));
		if (d_lp > chain->max_distance) return 0; 
		else poisson_dist_term (&sum_l, chain, chain->segment[t->prev->last_segment], &d_lp);
	}
	if (t->next) {
		d_rp = dSPR_topology_l (t->next, chain->topol->bkp[BKPshuffle_random], chain->split, &(chain->max_dist_split));
		if (d_rp > chain->max_distance) return 0; 
		else poisson_dist_term (&sum_r, chain, chain->segment[t->last_segment], &d_rp);
	}

	copy_topology_from_topology (chain->topol->bkp[BKPshuffle_copy], t);
	copy_topology_from_topology_mapping (t, chain->topol->bkp[BKPshuffle_random], chain->split);
	swap_likelihood (t, chain->segment, chain->topol->bkp[BKPshuffle_random]);
	
	chain->ln_likelihood_moved_branches (chain->segment, t);
//	if (debug_topol (t, chain)) printf ("update_topol_shuffle\n");

	acceptance = chain->kT * (t->likelihood_proposal - t->likelihood_current + sum_l + sum_r);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		accept_likelihood_moved_branches (chain->segment, t);

		if (t->prev) chain->segment[t->prev->last_segment]->dist = d_lp;
		if (t->next) chain->segment[t->last_segment      ]->dist = d_rp;  
		return 1;
	}
	else {
		unswap_likelihood (t, chain->segment, chain->topol->bkp[BKPshuffle_random]);
		copy_topology_from_topology (t, chain->topol->bkp[BKPshuffle_copy]);
		clear_topology (t);
		return 0;
	}
}

void
mcmc_topol_update_outside (chain_data chain, topology t, acceptance_frequency *afreq)
{
  double sum_[2], acceptance, transition = 0., 
         original_temp = chain->kT, temperature = (chain->kT * chain->kT_mini);
	int d_[4], i;

	sum_[0] = sum_[1] = 0.;
	d_[0] =	d_[1] =	d_[2] = d_[3] = 0;

	afreq->propose[FREQtopol]++;
	copy_topology_from_topology (chain->topol->bkp[BKPcopy], t);
	link_current_to_accepted (t, chain->segment);

	/* pi_beta(x)/pi(x) */
	if (t->prev) {
		d_[0] = chain->segment[t->prev->last_segment]->dist;
		poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->prev->last_segment], &d_[0], &(temperature), &(chain->kT));
	}
	if (t->next) {
		d_[1] = chain->segment[t->last_segment]->dist;
		poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->last_segment], &d_[1], &(temperature), &(chain->kT));
	}
	acceptance = ( ((temperature - chain->kT) * t->likelihood_current) + sum_[0] + sum_[1] );

	/* Al-wadhi-like mini-sampler (where x' and x are the same): heat/cold chain swap */
	chain->kT = temperature; 
	afreq->propose[FREQheat] += (chain->n_mini);
	for (i=0; i < chain->n_mini; i++) if ( mcmc_topol_update (chain, t)) {
		afreq->accept[FREQheat]++;
	}
	chain->kT = original_temp;

	/* pi(x_star)/pi_beta(x_star) */
	sum_[0] = sum_[1] = 0.;
	if (t->prev) {
		d_[2] = chain->segment[t->prev->last_segment]->dist;
		poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->prev->last_segment], &d_[2], &(chain->kT), &(temperature) );
	}
	if (t->next) {
		d_[3] = chain->segment[t->last_segment]->dist;
		poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->last_segment], &d_[3], &(chain->kT), &(temperature) );
	}
	acceptance += ( ((chain->kT - temperature) * t->likelihood_current) + sum_[0] + sum_[1] );

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		afreq->accept[FREQtopol]++;
		link_accepted_to_current (t, chain->segment);
	}
	else {
		if (t->prev) chain->segment[t->prev->last_segment]->dist = d_[0];
		if (t->next) chain->segment[t->last_segment]->dist = d_[1];
		copy_topology_from_topology (t, chain->topol->bkp[BKPcopy]);
		clear_topology (t);
	}
}

void
al_awadhi_minisampler (chain_data chain, topology t, acceptance_frequency *afreq)
{
	int i;
	double original_temperature = chain->kT;

	chain->kT *= chain->kT_cycles;

	afreq->propose[FREQcycle] += (chain->n_cycles);
	for (i=0; i < chain->n_cycles; i++) if ( mcmc_topol_update (chain, t)) {
		afreq->accept[FREQcycle]++;
	}

	chain->kT = original_temperature;
}

void
proposal_update_topol_birth_forward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2)
{
  double sum_[2], acceptance, transition, temperature = (chain->kT * chain->kT_cycles);
	int d_[4]; 
	
  afreq->propose[FREQadd]++;
	/* distances in order left --> right (where 0,1 are new dists and 2,3 are old */
	sum_[0] = sum_[1] = 0.;
	d_[0] =	d_[1] =	d_[2] = d_[3] = 0;

	/* this function will create t->prev and sort a random topology for it */
	add_topology_forward (t, chain->segment, chain->topol, chain->split);

  /* store new and old (if rejected in the end) dist info and prior ratio terms */
	if (t->prev->prev) {
		d_[0] = dSPR_topology_l (t->prev, t->prev->prev, chain->split, &(chain->max_dist_split));
		if (d_[0] > chain->max_distance) { reject_add_topology_forward (t, chain->segment, chain->topol); return; }
    poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->prev->prev->last_segment], &d_[0], &(temperature), &(chain->kT));
		d_[2] = chain->segment[t->prev->prev->last_segment]->dist;
  }

  d_[1] = dSPR_topology_l (t, t->prev, chain->split, &(chain->max_dist_split));
  if (d_[1] > chain->max_distance) { reject_add_topology_forward (t, chain->segment, chain->topol); return; }
  poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->prev->last_segment],  &d_[1], &(temperature), &(chain->kT));
	d_[3] = chain->segment[t->prev->last_segment]->dist;

  /* likelihood ratio of dimension transition */
	/* pi_beta(x')/pi(x) */
  chain->ln_likelihood_moved_branches (chain->segment, t->prev);
//	if (debug_topol (t->prev, chain)) printf ("birth_forward\n");
  acceptance = ( (temperature * t->prev->likelihood_proposal) - (chain->kT * t->prev->likelihood_current) );
	acceptance += (sum_[0] + sum_[1]);

	/* accept transient state */
	accept_likelihood_moved_branches (chain->segment, t->prev);
	if (t->prev->prev) chain->segment[t->prev->prev->last_segment]->dist = d_[0]; 
	chain->segment[t->prev->last_segment]->dist = d_[1];

	al_awadhi_minisampler (chain, t->prev, afreq);

	/* likelihood ratio of temperatures */
	/* pi(x_star)/pi_beta(x_star) */
	sum_[0] = sum_[1] = 0.;
  poisson_dist_term_no_ratio (&sum_[1], chain, chain->segment[t->prev->last_segment]);
  if (t->prev->prev) 
    poisson_dist_term_no_ratio (&sum_[0], chain, chain->segment[t->prev->prev->last_segment]);
  
	sum_[0] += (sum_[1] + t->prev->likelihood_current);
  acceptance += ( (chain->kT - temperature) * sum_[0] );
//	acceptance += ( chain->logD + log((double)(t->last_segment - t->prev->first_segment)) );
  transition = (*chain->logD) + (*log2);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		afreq->accept[FREQadd]++;
		accept_add_topology_forward (t, chain->segment);
		return;
	}
	else {
    chain->segment[t->prev->last_segment]->dist = d_[3];
    if (t->prev->prev) chain->segment[t->prev->prev->last_segment]->dist = d_[2];
		reject_add_topology_forward (t, chain->segment, chain->topol);
		return;
	}
}

void
proposal_update_topol_birth_backward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2)
{
  double sum_[2], acceptance, transition, temperature = (chain->kT * chain->kT_cycles);
	int d_[4]; 
	
  afreq->propose[FREQadd]++;
	/* distances in order right --> left (where 0,1 are new dists and 2,3 are old */
	sum_[0] = sum_[1] = 0.;
	d_[0] =	d_[1] =	d_[2] = d_[3] = 0;

	/* this function will create t->next and sort a random topology for it */
	add_topology_backward (t, chain->segment, chain->topol, chain->split);

	/* store new and old (if rejected in the end) dist info and prior ratio terms */
  if (t->next->next) {
    d_[0] = dSPR_topology_l (t->next->next, t->next, chain->split, &(chain->max_dist_split));
    if (d_[0] > chain->max_distance) { reject_add_topology_backward (t, chain->segment, chain->topol); return; }
    poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->next->last_segment], &d_[0], &(temperature), &(chain->kT));
    d_[2] = chain->segment[t->next->last_segment]->dist;
  }

  d_[1] = dSPR_topology_l (t->next, t, chain->split, &(chain->max_dist_split));
  if (d_[1] > chain->max_distance) { reject_add_topology_backward (t, chain->segment, chain->topol); return; }
  poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->last_segment],  &d_[1], &(temperature), &(chain->kT));
  d_[3] = chain->segment[t->last_segment]->dist;

  /* likelihood ratio of dimension transition */
	/* pi_beta(x')/pi(x) */
	chain->ln_likelihood_moved_branches (chain->segment, t->next);
//	if (debug_topol (t->next, chain)) printf ("birth_backward\n");
	acceptance  = ( (temperature * t->next->likelihood_proposal) - (chain->kT * t->next->likelihood_current) );
	acceptance += (sum_[0] + sum_[1]);

	/* accept transient state */
	accept_likelihood_moved_branches (chain->segment, t->next);
	if (t->next->next) chain->segment[t->next->last_segment]->dist = d_[0]; 
  chain->segment[t->last_segment]->dist = d_[1];

	al_awadhi_minisampler (chain, t->next, afreq);

	/* likelihood ratio of temperatures */
	/* pi(x_star)/pi_beta(x_star) */
	sum_[0] = sum_[1] = 0.;
	poisson_dist_term_no_ratio (&sum_[1], chain, chain->segment[t->last_segment]);
	if (t->next->next) 
		poisson_dist_term_no_ratio (&sum_[0], chain, chain->segment[t->next->last_segment]);

	sum_[0] += (sum_[1] + t->next->likelihood_current);
	acceptance += ( (chain->kT - temperature) * sum_[0] );
//	acceptance += ( chain->logD + log((double)(t->next->last_segment - t->first_segment)) );
  transition = (*chain->logD) + (*log2);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		afreq->accept[FREQadd]++;
		accept_add_topology_backward (t, chain->segment);
    return;
  }
  else {
		chain->segment[t->last_segment]->dist = d_[3];
		if (t->next->next) chain->segment[t->next->last_segment]->dist = d_[2];
		reject_add_topology_backward (t, chain->segment, chain->topol);
    return;
  }
}

void
proposal_update_topol_death_forward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2)
{
  double sum_[2], acceptance, transition, temperature = (chain->kT * chain->kT_cycles);
  int d_[4]; 
	
	afreq->propose[FREQremove]++;
	/* distances in order left --> right where (0,1) are OLD and (2,3) are NEW dists */
	sum_[0] = sum_[1] = 0.;
	d_[0] =	d_[1] =	d_[2] = d_[3] = 0;

	/* likelihood ratio of temperatures and store old distance info */
	/* pi_beta(x_star)/pi(x_star) */
	if (t->next) { 
		d_[3] = dSPR_topology_l (t->next, t->prev, chain->split, &(chain->max_dist_split));
		if (d_[3] > chain->max_distance) { return; }
		poisson_dist_term_no_ratio (&sum_[1], chain, chain->segment[t->last_segment]);
		d_[1] = chain->segment[t->last_segment]->dist;
	}
	poisson_dist_term_no_ratio (&sum_[0], chain, chain->segment[t->prev->last_segment]);
	d_[0] = chain->segment[t->prev->last_segment]->dist;
	
	sum_[0] += (sum_[1] + t->likelihood_accepted);
	acceptance = ( (temperature - chain->kT) * sum_[0] );

	pre_add_remove_topology (t, chain->segment, chain->topol);
	al_awadhi_minisampler (chain, t, afreq);
	remove_topology_forward (t, chain->segment, chain->topol, chain->split);

	/* new dist info (t and t->prev are the same) and prior ratio */
	sum_[0] = sum_[1] = 0.;
	poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->prev->last_segment], &d_[2], &(chain->kT), &(temperature));
	if (t->next) {
		poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->last_segment], &d_[3], &(chain->kT), &(temperature));
	}

	/* likelihood ratio of dimension transition */
	/* pi(x)/pi_beta(x') */
	chain->ln_likelihood_moved_branches (chain->segment, t);
//	if (debug_topol (t, chain)) printf ("death_forward\n");
	acceptance += ( (chain->kT * t->likelihood_proposal) - (temperature * t->likelihood_current) );
	acceptance += (sum_[0] + sum_[1]);
//	acceptance -= ( chain->logD + log((double)(t->last_segment - t->prev->first_segment)) );
  transition = (*log2) - (*chain->logD);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		if (t->next) chain->segment[t->last_segment]->dist = d_[3];
		chain->segment[t->prev->last_segment]->dist = d_[2];
		accept_likelihood_moved_branches (chain->segment, t);
    afreq->accept[FREQremove]++;
		accept_remove_topology_forward (t, chain->segment, chain->topol);
		return;
	}
	else {
		if (t->next) chain->segment[t->last_segment]->dist = d_[1];
		chain->segment[t->prev->last_segment]->dist = d_[0];
		reject_remove_topology_forward (t, chain->segment, chain->topol);
		return;
	}
}

void
proposal_update_topol_death_backward (topology t, chain_data chain, acceptance_frequency *afreq, double *log2)
{
  double sum_[2], acceptance, transition, temperature = (chain->kT * chain->kT_cycles);
	int d_[4]; 
	
	afreq->propose[FREQremove]++;
	/* distances in order right --> left where (0,1) are OLD and (2,3) are NEW dists */
	sum_[0] = sum_[1] = 0.;
	d_[0] =	d_[1] =	d_[2] = d_[3] = 0;


  /* likelihood ratio of temperatures and store old distance info */
	/* pi_beta(x_star)/pi(x_star) */
	if (t->prev) { 
		d_[3] = dSPR_topology_l (t->next, t->prev, chain->split, &(chain->max_dist_split));
    if (d_[3] > chain->max_distance) { return; }
		poisson_dist_term_no_ratio (&sum_[1], chain, chain->segment[t->prev->last_segment]);
		d_[1] = chain->segment[t->prev->last_segment]->dist;
	}
	poisson_dist_term_no_ratio (&sum_[0], chain, chain->segment[t->last_segment]);
	d_[0] = chain->segment[t->last_segment]->dist;
	sum_[0] += (sum_[1] + t->likelihood_accepted);
	acceptance = ( (temperature - chain->kT) * sum_[0] );

	pre_add_remove_topology (t, chain->segment, chain->topol);
	al_awadhi_minisampler (chain, t, afreq);
	remove_topology_backward (t, chain->segment, chain->topol, chain->split);

  /* new dist info (t and t->next are the same) and prior ratio */
	sum_[0] = sum_[1] = 0.;
	poisson_dist_term_awadhi (&sum_[0], chain, chain->segment[t->last_segment], &d_[2], &(chain->kT), &(temperature));
	if (t->prev) {
		poisson_dist_term_awadhi (&sum_[1], chain, chain->segment[t->prev->last_segment], &d_[3], &(chain->kT), &(temperature));
	}

	/* likelihood ratio of dimension transition */
	/* pi(x)/pi_beta(x') */
	chain->ln_likelihood_moved_branches (chain->segment, t);
//	if (debug_topol (t, chain)) printf ("death_backward\n");
	acceptance += ( (chain->kT * t->likelihood_proposal) - (temperature * t->likelihood_current) );
	acceptance += (sum_[0] + sum_[1]);
//	acceptance -= ( chain->logD + log((double)(t->next->last_segment - t->first_segment)) );
  transition = (*log2) - (*chain->logD);

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		if (t->prev) chain->segment[t->prev->last_segment]->dist = d_[3];
		chain->segment[t->last_segment]->dist = d_[2];
		accept_likelihood_moved_branches (chain->segment, t);
    afreq->accept[FREQremove]++;
		accept_remove_topology_backward (t, chain->segment, chain->topol);
		return;
  }
  else {
    if (t->prev) chain->segment[t->prev->last_segment]->dist = d_[1];
    chain->segment[t->last_segment]->dist = d_[0];
    reject_remove_topology_backward (t, chain->segment, chain->topol);
		return;
	}
}

void
proposal_update_topol_shift_forward (topology t, chain_data chain, acceptance_frequency *afreq)
{
  double sum_l = 0., sum_r = 0., acceptance, transition;
	int d_lp = 0, d_rp = chain->segment[ t->prev->last_segment ]->dist;
	afreq->propose[FREQshift]++;
	
	shift_topology_forward (t, chain->segment, chain->topol, chain->split);

	poisson_dist_term (&sum_l, chain, chain->segment[t->prev->prev->last_segment], &d_lp);
	poisson_dist_term (&sum_r, chain, chain->segment[t->prev->last_segment      ], &d_rp);
	chain->ln_likelihood_moved_branches (chain->segment, t->prev);
//	if (debug_topol (t->prev, chain)) printf ("shift_forward\n");
	
	acceptance = chain->kT * (t->prev->likelihood_proposal - t->prev->likelihood_current + sum_l + sum_r);
  transition = log((double)(t->last_segment - t->prev->first_segment)); 
  transition -= log((double)(t->prev->last_segment - t->prev->prev->first_segment));

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		accept_likelihood_moved_branches (chain->segment, t->prev);
		chain->segment[t->prev->prev->last_segment]->dist = d_lp;
		chain->segment[t->prev->last_segment      ]->dist = d_rp;  
		afreq->accept[FREQshift]++;
		accept_shift_topology_forward (t, chain->segment, chain->topol);
	}
	else {
		reject_shift_topology_forward (t, chain->segment, chain->topol);
	}
}

void
proposal_update_topol_shift_backward (topology t, chain_data chain, acceptance_frequency *afreq)
{	
  double sum_l = 0., sum_r = 0., acceptance, transition;
	int d_lp = chain->segment[ t->last_segment ]->dist, d_rp = 0; 
	afreq->propose[FREQshift]++;

	shift_topology_backward (t, chain->segment, chain->topol, chain->split);

	poisson_dist_term (&sum_l, chain, chain->segment[t->last_segment      ], &d_lp);
	poisson_dist_term (&sum_r, chain, chain->segment[t->next->last_segment], &d_rp);
	chain->ln_likelihood_moved_branches (chain->segment, t->next);
//	if (debug_topol (t->next, chain)) printf ("shift_forward\n");

	acceptance = chain->kT * (t->next->likelihood_proposal - t->next->likelihood_current + sum_l + sum_r);
  transition = log((double)(t->next->last_segment - t->first_segment)); 
  transition -= log((double)(t->next->next->last_segment - t->next->first_segment));

  if (chain->acceptance_probability (chain, &acceptance, &transition)) {
		accept_likelihood_moved_branches (chain->segment, t->next);
		chain->segment[t->last_segment      ]->dist = d_lp;
		chain->segment[t->next->last_segment]->dist = d_rp;
		afreq->accept[FREQshift]++;
		accept_shift_topology_backward (t, chain->segment, chain->topol);
	}
	else {
		reject_shift_topology_backward (t, chain->segment, chain->topol);
	}
}

bool
debug_topol (topology t, chain_data chain)
{
	double old_lik = t->likelihood_proposal;
	int i;
	/*check likelihood */
	ln_likelihood (chain->segment, t);

	if ((t->likelihood_proposal - old_lik)*(t->likelihood_proposal - old_lik) > 0.1) {
		for (i=t->nleaves; i < t->nnodes; i++) {
			printf ("%02d  %02d  %02d ", t->nodelist[i]->id, t->nodelist[i]->map_id, t->nodelist[i]->up->id);
			printf (" [ %d %d ]\n", (int) t->nodelist[i]->u_done, (int) t->nodelist[i]->d_done);
		}
		printf ("DDBUG ( %6.2f %6.2f ) \t %d ", t->likelihood_proposal, old_lik, t->id);
		printf (". %d %d . ", t->first_segment, t->last_segment);
		return true;
	}
	return false;
}

