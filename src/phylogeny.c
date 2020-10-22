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

/*! \file phylogeny.c 
 *  \brief 
 */

#include "phylogeny.h"

node_likelihood new_node_likelihood (int nchar, int ncycles);
void store_likelihood_at_leaf (node_likelihood l, char *align, int *index, int align_size);

phylogeny
new_phylogeny (int nchar, int ntax, int ncycles)
{
	phylogeny p;
	int i;

	p = 
	(phylogeny) biomc2_malloc (sizeof (struct phylogeny_struct));
	p->nchar  = nchar;
	p->ntax   = ntax;
	p->nnodes = ((2*p->ntax) - 2);
	p->likelihood_accepted = p->likelihood_current = p->likelihood_proposal = 0.;
	p->rate = -1.; /* defined by tree or sequence in chain_data */
	p->dist = 0;
	p->kappa = p->lambda = p->penalty = 0.;

	p->l = 
	(node_likelihood*) biomc2_malloc (p->nnodes * sizeof (node_likelihood));
	
	p->l_ptr = 
	(node_likelihood*) biomc2_malloc ((p->nnodes - p->ntax) * sizeof (node_likelihood));
	
	/* allocate memory for likelihood vectors */
	for (i=0; i < p->ntax; i++)         p->l[i] = new_node_likelihood (p->nchar, 1); 
	for (i=p->ntax; i < p->nnodes; i++) p->l[i] = new_node_likelihood (p->nchar, ncycles + 3); 

	for (i=p->ntax; i < p->nnodes; i++) p->l_ptr[i - p->ntax] = p->l[i];

	p->site_weight = (double*) biomc2_malloc (p->nchar * sizeof (double));
	for (i=0; i < p->nchar; i++) p->site_weight[i] = 0.; 

  return p;
}

node_likelihood
new_node_likelihood (int nchar, int ncycles)
{
	int i,j;
	node_likelihood l;
	l =	
	(node_likelihood) biomc2_malloc (sizeof (struct node_likelihood_struct));
	/* one vector for current state, one for proposal update and one for accepted */
	l->size = ncycles; 
	l->u = 
	(lk_vector*) biomc2_malloc ((l->size) * sizeof (lk_vector));
	l->d = 
	(lk_vector*) biomc2_malloc ((l->size) * sizeof (lk_vector));
	for (i=0; i < l->size; i++) {
		l->u[i] = 
		(lk_vector) biomc2_malloc (sizeof (struct lk_vector_struct));
		l->d[i] = 
		(lk_vector) biomc2_malloc (sizeof (struct lk_vector_struct));
	}

	for (i=0; i < l->size; i++) {
		l->u[i]->lk = 
		(double**) biomc2_malloc ((nchar) * sizeof (double*));
		l->d[i]->lk = 
		(double**) biomc2_malloc ((nchar) * sizeof (double*));
		for (j=0; j < (nchar); j++) {
			l->u[i]->lk[j] =
			(double*) biomc2_malloc (4 * sizeof (double));
			l->d[i]->lk[j] =
			(double*) biomc2_malloc (4 * sizeof (double));
		}
		l->u[i]->prev = l->u[i]->next = 
		l->d[i]->prev = l->d[i]->next = NULL;
	}
	
	l->u[0]->prev = l->u[l->size-1];
	l->d[0]->prev = l->d[l->size-1];
	l->u[l->size-1]->next = l->u[0];
	l->d[l->size-1]->next = l->d[0];

	for (i=0; i < l->size; i++) {
		if (!l->u[i]->prev) l->u[i]->prev = l->u[i-1];
		if (!l->d[i]->prev) l->d[i]->prev = l->d[i-1];
		if (!l->u[i]->next) l->u[i]->next = l->u[i+1];
		if (!l->d[i]->next) l->d[i]->next = l->d[i+1];
	}

	l->u_current = l->u_accepted = l->u[0];
	l->d_current = l->d_accepted = l->d[0];

	return l;
}

void 
del_phylogeny (phylogeny p)
{
	int i, j, k, nnodes = 2*p->ntax - 2;
	
	if (!p) return;

	if (p->l) { 
		for (i=nnodes-1; i >=0; i--) {
			if (p->l[i] && p->l[i]->u && p->l[i]->d) {
				for (k= (p->l[i]->size - 1); k>=0; k--) {
					if (p->l[i]->u[k] && p->l[i]->u[k]->lk) 
					 {
						for (j=p->nchar-1; j>=0; j--) if (p->l[i]->u[k]->lk[j]) free (p->l[i]->u[k]->lk[j]);
						free (p->l[i]->u[k]->lk); 
						free (p->l[i]->u[k]); 
					 }
					if (p->l[i]->d[k] && p->l[i]->d[k]->lk) 
					 { 
						for (j=p->nchar-1; j>=0; j--) if (p->l[i]->d[k]->lk[j]) free (p->l[i]->d[k]->lk[j]);
						free (p->l[i]->d[k]->lk);
						free (p->l[i]->d[k]);
					 }
				}
				free (p->l[i]->u);
				free (p->l[i]->d);
				free (p->l[i]);
			}
		}
		free (p->l);
		free (p->l_ptr);
	}
	if (p->site_weight) free (p->site_weight);
	free (p);
}

phylogeny
new_phylogeny_from_nexus (nexus_alignment a, double *siteweight, int *idx, int n_idx, int ncycles)
{
  phylogeny p;
  int i;

	p = new_phylogeny (n_idx, a->ntax, ncycles);
	for (i=0; i < n_idx; i++) p->site_weight[i] = siteweight[ idx[i] ];
	for (i=0; i < a->ntax; i++) store_likelihood_at_leaf (p->l[i], a->character[i], idx, n_idx);
  return p;
}

void
store_likelihood_at_leaf (node_likelihood l, char *align, int *index, int align_size)
{
	int i, j;
	int bit[256];

	for (j=0; j < align_size; j++)
		for (i=0; i < 4; i++)
			l->u[0]->lk[j][i] = l->d[0]->lk[j][i] = 0.;

  for (i = 0; i < 256; i++) bit[i] = 0;

  bit['A'] = 1;    /* .   A */
  bit['B'] = 14;   /* .TGC  */
  bit['C'] = 2;    /* .  C  */
  bit['D'] = 13;   /* .TG A */
  bit['G'] = 4;    /* . G   */
  bit['H'] = 11;   /* .T CA */
  bit['K'] = 12;   /* .TG   */
  bit['M'] = 3;    /* .  CA */
  bit['N'] = 15;   /* .TGCA */
  bit['O'] = 15;   /* .TGCA */
  bit['R'] = 5;    /* . G A */
  bit['S'] = 6;    /* . GC  */
  bit['T'] = 8;    /* .T    */
  bit['U'] = 8;    /* .T    */
  bit['V'] = 7;    /* . GCA */
  bit['W'] = 9;    /* .T  A */
  bit['X'] = 15;   /* .TGCA */
  bit['Y'] = 10;   /* .T C  */
  bit['?'] = 15;   /* .TGCA */
  bit['-'] = 15;   /* .TGCA */
//bit['-'] = 16;   /* -TGCA */
  
	for (j = 0; j < align_size; j++)
		for (i = 0; i < 4; i++)
			if (bit[(int)(align[ index[j] ])] & (1 << i)) l->u[0]->lk[j][i] = l->d[0]->lk[j][i] = 1.;

}

double
average_rate_distance_method (phylogeny phy)
{
	int i1, i2, b1, b2, pos;
	int e_1=0, d_1=0;
	double e_2 = 0., d_2=0.;
	double **seq1, **seq2;
	
	for (i1=1; i1 < phy->ntax; i1++) {
		seq1 = phy->l[i1]->d[0]->lk;
		for (i2=0; i2 < i1; i2++) {
			seq2 = phy->l[i2]->d[0]->lk;

			for (pos=0; pos < phy->nchar; pos++) {
				e_1 = d_1 = 0;

        for (b1=0; b1 < 4; b1++) for (b2=0; b2 < 4; b2++) if ((seq1[pos][b1] * seq2[pos][b2]) > 1e-6) {
          if (b1 == b2) e_1++;
          else d_1++;
        }
        e_2 += phy->site_weight[pos] * (double)(e_1);
        d_2 += phy->site_weight[pos] * (double)(d_1);
      }
	
		}
	}
      
  e_2 += d_2; /* total number of observations */
  d_2 /= e_2; /* fraction of differences */
	d_2 /= (phy->ntax);
	
	if (d_2 < 0.75) { /* Jukes-Cantor distance */
		d_2 *= 4./3.; d_2 = log (1. - d_2); d_2 *= -3./4.;
	}
  else d_2 = 1;

	return d_2;
}

void
update_equilibrium_frequencies (double *pi, phylogeny phy)
{
	int i, b, pos;
	double **seq;
	
	for (i=0; i < phy->ntax; i++) {
		seq = phy->l[i]->d[0]->lk;
		for (pos=0; pos < phy->nchar; pos++) 
			for (b=0; b < 4; b++)
				pi[b] += (phy->site_weight[pos] * seq[pos][b]);
	}
}

