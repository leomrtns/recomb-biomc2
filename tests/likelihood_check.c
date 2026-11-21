/* Compare full, cached and remapped likelihoods with independent long-double pruning. */
#include "chain_data.h"

static long double
log_transition (const double q[4], const long double state[4])
{
  long double largest = state[0], sum = 0.;
  for (int i=1; i < 4; i++) if (state[i] > largest) largest = state[i];
  if (!isfinite (largest)) return largest;
  for (int i=0; i < 4; i++) sum += q[i] * expl (state[i] - largest);
  return logl (sum) + largest;
}

static void
reference_subtree (topol_node node, phylogeny p, int site, long double result[4])
{
  if (!node->internal) {
    for (int i=0; i < 4; i++) result[i] = logl (p->l[node->id]->d_accepted->lk[site][i]);
    return;
  }
  long double left[4], right[4];
  reference_subtree (node->left, p, site, left);
  reference_subtree (node->right, p, site, right);
  for (int i=0; i < 4; i++) {
    result[i] = log_transition (p->Q[i], left) + log_transition (p->Q[i], right);
  }
}

static double
reference_likelihood (topology t, phylogeny p)
{
  long double result = 0.;
  for (int site=0; site < p->nchar; site++) {
    long double lower[4], likelihood = 0.;
    reference_subtree (t->root->left, p, site, lower);
    long double root[4], largest;
    for (int i=0; i < 4; i++)
      root[i] = logl (p->pi[i] * p->l[0]->u_accepted->lk[site][i]) + log_transition (p->Q[i], lower);
    largest = root[0];
    for (int i=1; i < 4; i++) if (root[i] > largest) largest = root[i];
    for (int i=0; i < 4; i++) likelihood += expl (root[i] - largest);
    result += p->site_weight[site] * (logl (likelihood) + largest);
  }
  return (double)result;
}

static void
check_score (double score, double expected)
{
  if (!isfinite (score) || fabs (score - expected) > 1.e-8 * (1. + fabs (expected)))
    biomc2_error ("likelihood mismatch: %.15g versus %.15g", score, expected);
}

int
main (int argc, char **argv)
{
  if (argc != 2) return EXIT_FAILURE;
  random_number = new_biomc2_rng ();
  random_number->s1 = 1234567; random_number->s2 = 2345678; random_number->s3 = 3456789;
  nexus_treespace trees = read_nexus_treespace_file (argv[1]);
  topology t = new_topology (trees->nleaves), old = new_topology (trees->nleaves);
  topology proposal = new_topology (trees->nleaves);
  split_space split = new_split_space (trees->nleaves);
  copy_topology_from_nexus_tree (t, trees->T[0]);
  t->first_segment = t->last_segment = 0;
  double pi[] = {.1, .2, .3, .4};
  phylogeny p = new_phylogeny (4, trees->nleaves, 2);
  p->pi = pi;
  for (int i=0; i < 4; i++) {
    p->site_weight[i] = i + 1;
    for (int j=0; j < 4; j++) p->Q[i][j] = .85 * (i == j) + .15 * pi[j];
  }
  for (int taxon=0; taxon < trees->nleaves; taxon++) for (int site=0; site < 4; site++) {
    for (int state=0; state < 4; state++)
      p->l[taxon]->u_accepted->lk[site][state] = p->l[taxon]->d_accepted->lk[site][state] =
        site == 3 || state == (taxon + site) % 4;
    p->l[taxon]->u_accepted->lk[site][4] = p->l[taxon]->d_accepted->lk[site][4] = 0.;
  }
  ln_likelihood_proposal (p, t);
  check_score (p->likelihood_proposal, reference_likelihood (t, p));
  accept_likelihood_proposal (p, t);
  link_current_to_accepted (t, &p);
  for (int move=0; move < 40; move++) {
    if (move == 20) {
      for (int i=0; i < 4; i++) for (int j=0; j < 4; j++)
        p->Q[i][j] = .2 * (i == j) + .8 * pi[j];
      ln_likelihood (&p, t);
      check_score (t->likelihood_proposal, reference_likelihood (t, p));
      ln_likelihood_proposal (p, t);
      check_score (p->likelihood_proposal, reference_likelihood (t, p));
      accept_likelihood_proposal (p, t);
      link_current_to_accepted (t, &p);
    }
    copy_topology_from_topology (old, t);
    bool remap = move % 5 == 0;
    if (remap) {
      copy_topology_from_topology (proposal, t);
      for (int i=0; i < 3; i++) spr_apply_random_move_spr (proposal, true);
      copy_topology_from_topology_mapping (t, proposal, split);
      swap_likelihood (t, &p, proposal);
    }
    else if (move % 2) spr_apply_random_move_spr (t, false);
    else spr_apply_random_move_nni (t, false);
    double expected = reference_likelihood (t, p);
    /* A full recalculation must work before a cached proposal has filled next. */
    ln_likelihood (&p, t);
    check_score (t->likelihood_proposal, expected);
    ln_likelihood_moved_branches (&p, t);
    check_score (t->likelihood_proposal, expected);
    if (move % 3) {
      accept_likelihood_moved_branches (&p, t);
      link_accepted_to_current (t, &p);
    }
    else {
      if (remap) unswap_likelihood (t, &p, proposal);
      copy_topology_from_topology (t, old);
    }
    clear_topology (t);
  }
  del_phylogeny (p);
  del_topology (t); del_topology (old); del_topology (proposal);
  del_split_space (split);
  del_nexus_treespace (trees);
  del_biomc2_rng (random_number);
  return EXIT_SUCCESS;
}
