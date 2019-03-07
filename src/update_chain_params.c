#include "update_chain.h"

/* local function prototypes */

void update_variable_logscale (double *result, double *current, double *step_size);
void update_variable (double *result, double *current, double *step_size);

void
full_likelihood (double *r, chain_data chain)
{
  int i;
	double d1 = 0.;
	topology t;

	(*r) = ( (chain->prior_rate/chain->hp_rate) * log (chain->hp_rate) );
	for (t=chain->topol->first; (t); t = t->next) (*r) += t->likelihood_accepted;
	for (i=0; i < chain->n_segments; i++) d1 += chain->segment[i]->rate;
	(*r) += ( (d1/chain->prior_rate) * log (chain->prior_rate) );

  for (i=0; i < chain->n_segments - 1; i++) {
		d1 = chain->segment[i]->penalty + WEIGHT;
		(*r) -= ( chain->lnN[i] + d1 * chain->factorial[ chain->segment[i]->dist ] ) ;
		(*r) += ( log (chain->segment[i]->lambda) * (chain->hp_lambda_alpha - 1. + (chain->segment[i]->dist * d1)) );
		(*r) += ( (chain->hp_penalty_alpha - 1.) * log (chain->segment[i]->penalty) );
		(*r) -= ( chain->segment[i]->lambda * (d1 + chain->hp_lambda_beta) + 
							chain->hp_penalty_beta * chain->segment[i]->penalty );
	}
}

void 
lnTruncate (double *r, double *w, double *l, chain_data chain)
{
  int j;
  double lg, lnl = log (*l);

  (*r) = 0.;
  for (j=0; j <= chain->max_distance; j++) {
    lg = (double)(j) * lnl - (*l) -  chain->factorial[j];
    (*r) += exp ((*w) * lg);
  }
  (*r) = log ((*r));
}

void 
calculate_moments (double *r1, double *r2, int i, chain_data chain)
{
	int j;
	double lg, lnl = log (chain->segment[i]->lambda), w = chain->segment[i]->penalty + WEIGHT;

	(*r1) = (*r2) = 0.;
	/* first moment about zero (E[X]) */
	for (j=1; j <= chain->max_distance; j++) {
		lg = (double)(j) * lnl - (chain->segment[i]->lambda) - chain->factorial[j];
		lg = exp (w * lg);
		(*r1) += ( (double)(j)     * lg ); 
	}
	(*r1) /= exp (chain->lnN[i]);
	/* second moment about E[X], also called second central moment (Var[X]) */
	for (j=1; j <= chain->max_distance; j++) {
		lg = (double)(j) * lnl - (chain->segment[i]->lambda) - chain->factorial[j];
		lg = exp (w * lg);
		(*r2) += ( ((double)(j)-(*r1)) * ((double)(j)-(*r1)) * lg ); 
	}
	(*r2) /= exp (chain->lnN[i]);
}

void
propose_swap_chains (chain_data chain1, chain_data chain2, acceptance_frequency *cold_freq, acceptance_frequency *heat_freq)
{
	double lnl_1, lnl_2;

	cold_freq->propose[FREQswapchain]++;
	heat_freq->propose[FREQswapchain]++;

	full_likelihood (&lnl_1, chain1);
	full_likelihood (&lnl_2, chain2);

	if (log (biomc2_rng_uniform_pos (random_number)) < ((chain2->kT - chain1->kT) * (lnl_1 - lnl_2)) ) {
		lnl_1      = chain1->kT;
		chain1->kT = chain2->kT;
		chain2->kT = lnl_1;

		cold_freq->accept[FREQswapchain]++;
		heat_freq->accept[FREQswapchain]++;
	}
}

void
update_variable_logscale (double *result, double *current, double *step_size)
{
	double u = biomc2_rng_uniform (random_number);
  (*result) = (*current) + (*step_size) * (u - 0.5);
}

void
update_variable (double *result, double *current, double *step_size)
{
  double u = biomc2_rng_uniform (random_number);
  (*result) = (*step_size) * (u - 0.5);
  (*result) = (*current) * exp (*result);
}

void
proposal_update_lambda (chain_data chain, int n, acceptance_frequency *afreq)
{
  double nlambda, ln_r, w = chain->segment[n]->penalty + WEIGHT, transition;
	
	afreq->propose[FREQlambda]++;

	update_variable (&nlambda, &(chain->segment[n]->lambda), &(chain->update_lambda));
	lnTruncate (&(chain->lnN_proposal[n]), &w, &nlambda, chain);

	ln_r  = (chain->segment[n]->lambda - nlambda) * (w + chain->hp_lambda_beta);
	ln_r += ( (w * (double)(chain->segment[n]->dist) + chain->hp_lambda_alpha - 1.) * 
						(log (nlambda) - log (chain->segment[n]->lambda)) );
	ln_r += ( chain->lnN[n] - chain->lnN_proposal[n] );

  ln_r = (chain->kT * ln_r);
  transition = log (nlambda) - log (chain->segment[n]->lambda);

  if (chain->acceptance_probability (chain, &ln_r, &transition)) {
		chain->lnN[n] = chain->lnN_proposal[n];
		chain->segment[n]->lambda = nlambda;
		afreq->accept[FREQlambda]++;
	}
}

void
proposal_update_penalty (chain_data chain, int n, acceptance_frequency *afreq)
{
  double npenalty, nw, r2 = 0., transition;
	
	afreq->propose[FREQpenalty]++;

	update_variable (&npenalty, &(chain->segment[n]->penalty), &(chain->update_penalty));
	nw = npenalty + WEIGHT;

	lnTruncate (&(chain->lnN_proposal[n]), &nw, &(chain->segment[n]->lambda), chain);
	r2 = ( chain->hp_penalty_beta + chain->factorial[ chain->segment[n]->dist ] + 
				 chain->segment[n]->lambda - ((double)(chain->segment[n]->dist) * log (chain->segment[n]->lambda)) );
	r2 *= ( chain->segment[n]->penalty - npenalty );
	r2 += ( chain->lnN[n] - chain->lnN_proposal[n] );
	r2 += ( (log (npenalty) - log (chain->segment[n]->penalty) ) * (chain->hp_penalty_alpha - 1.) );

	r2 = (chain->kT * r2);
  transition = log (npenalty) - log (chain->segment[n]->penalty);

  if (chain->acceptance_probability (chain, &r2, &transition)) {
		chain->segment[n]->penalty = npenalty;
		chain->lnN[n] = chain->lnN_proposal[n];
		afreq->accept[FREQpenalty]++;
	}
}

void 
proposal_update_substitution_rate (chain_data chain, int n, topology t, acceptance_frequency *afreq)
{
	double new_rate, a2, diff_likelihood, transition;

	afreq->propose[FREQrate]++;

  update_variable (&new_rate, &(chain->segment[n]->rate), &(chain->update_rate));
  update_Q_matrix (chain->segment[n], new_rate);

	chain->ln_likelihood_proposal (chain->segment[n], t);

	diff_likelihood = chain->segment[n]->likelihood_proposal - chain->segment[n]->likelihood_accepted;

  a2 = chain->kT * (diff_likelihood + (chain->segment[n]->rate - new_rate)/chain->prior_rate);

  transition = log (new_rate) - log (chain->segment[n]->rate);

  if (chain->acceptance_probability (chain, &a2, &transition)) {
		accept_likelihood_proposal (chain->segment[n], t);
		chain->segment[n]->rate = new_rate;
		t->likelihood_accepted += diff_likelihood;
		afreq->accept[FREQrate]++;
  }
	else update_Q_matrix (chain->segment[n], chain->segment[n]->rate);
}

void 
proposal_update_substitution_kappa (chain_data chain, int n, topology t, acceptance_frequency *afreq)
{
  double new_kappa, a2, diff_likelihood, transition;

	afreq->propose[FREQkappa]++;

	update_variable (&new_kappa, &(chain->segment[n]->kappa), &(chain->update_kappa));
	update_Q_eigenvalues (chain->segment[n], new_kappa);
	update_Q_matrix (chain->segment[n], chain->segment[n]->rate);

	chain->ln_likelihood_proposal (chain->segment[n], t);
	diff_likelihood = chain->segment[n]->likelihood_proposal - chain->segment[n]->likelihood_accepted;

  a2 = chain->kT * (diff_likelihood + (chain->segment[n]->kappa - new_kappa )/chain->prior_kappa);

  transition = log (new_kappa) - log (chain->segment[n]->kappa);

  if (chain->acceptance_probability (chain, &a2, &transition)) {
		accept_likelihood_proposal (chain->segment[n], t);
		chain->segment[n]->kappa = new_kappa;
		t->likelihood_accepted += diff_likelihood;
		afreq->accept[FREQkappa]++;
	}
	else {
		update_Q_eigenvalues (chain->segment[n], chain->segment[n]->kappa);
		update_Q_matrix (chain->segment[n], chain->segment[n]->rate);
	}
}

void
proposal_update_prior_rate (chain_data chain, acceptance_frequency *afreq)
{
  double new_prior, diff, sum = 0., transition;
  int i;
	
	afreq->propose[FREQrateprior]++;

	update_variable (&new_prior, &(chain->prior_rate), &(chain->update_rate));

  for (i=0; i < chain->n_segments; i++) sum += chain->segment[i]->rate;

	diff  = sum * (1./chain->prior_rate - 1./new_prior);
	diff += ( 1./chain->hp_rate * (chain->prior_rate - new_prior) );
	diff += ( (double)(chain->n_segments) * (log (chain->prior_rate) - log (new_prior)) );

  diff = (chain->kT * diff);
  transition = log (new_prior) - log (chain->prior_rate);

  if (chain->acceptance_probability (chain, &diff, &transition)) {
		chain->prior_rate = new_prior;
		afreq->accept[FREQrateprior]++;
	}
}

void
proposal_update_prior_kappa (chain_data chain, acceptance_frequency *afreq)
{
  double new_prior, diff, sum = 0., transition; 
	int i;

	afreq->propose[FREQkappaprior]++;

	update_variable (&new_prior, &(chain->prior_kappa), &(chain->update_kappa));

	for (i=0; i < chain->n_segments; i++) sum += chain->segment[i]->kappa;

	diff  = sum * (1./chain->prior_kappa - 1./new_prior);
	diff += ( 1./chain->hp_kappa * (chain->prior_kappa - new_prior) );
	diff += ( (double)(chain->n_segments) * (log (chain->prior_kappa) - log (new_prior)) );

  diff = (chain->kT * diff);
  transition = log (new_prior) - log (chain->prior_kappa);

  if (chain->acceptance_probability (chain, &diff, &transition)) {
		chain->prior_kappa = new_prior;
		afreq->accept[FREQkappaprior]++;
	}
}

