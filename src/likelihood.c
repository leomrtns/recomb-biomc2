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

/*! \file likelihood.c 
 *  \brief likelihod calculation functions
 */

#include "likelihood.h"

/* local variables */

int s1, s2, site;
double lnLikSequence = 0., LikSite,
			 lkl, lkr,
			 **up, **down, **sister, **left, **right;

/* local function prototypes */
void          update_likelihood (phylogeny *p, topol_node this, topology t);
void update_accepted_likelihood (phylogeny *p, topol_node this, topol_node upnode, topology t);
void          update_likelihood_moved_branches (phylogeny *p, topol_node this, topology t);
void update_accepted_likelihood_moved_branches (phylogeny *p, topol_node this, topol_node upnode, topology t);
void          update_likelihood_proposal (phylogeny p, topol_node this);
void update_accepted_likelihood_proposal (phylogeny p, topol_node this, topol_node upnode);
void mark_d_done (phylogeny *p, topol_node this, topology t);
void mark_d_done_moved_branches (phylogeny *p, topol_node this, topology t);


void
update_Q_eigenvalues (phylogeny p, double kappa)
{
  double k = ( kappa * ((p->pi[0]*p->pi[2]) + (p->pi[1]*p->pi[3])) + (p->pi[4]*p->pi[5]) );
//  k = (2. * k)/(2. + kappa);
  k *= 2.;
	p->psi[0] = 0.;
	p->psi[1] = 1./k;
	p->psi[2] = ((kappa * p->pi[5]) + p->pi[4])/k;
	p->psi[3] = ((kappa * p->pi[4]) + p->pi[5])/k;
}

/*! \brief  HKY model integrated over branch length 
 *
 *  Calculates \f$ Prob(j/i,\lambda) \f$ where \f$ i,j\in {A,C,G,T}\f$ and branch length 
 *  is assumed to have exponential distribution with mean \f${\lambda}\f$, leading
 *  to
 *  \f[ Prob(j/i,\lambda)=\int Prob(j/i,t)p(t/\lambda)\, \mathrm{d}t \f]
 *  Assuming 
 *  \f[ p(t/\lambda)\mathrm{d}t=\frac1{\lambda} e^{\frac{t}{\lambda}}\mathrm{d}t \f]
 *  and the HKY model \f$ Prob(j/i)\f$ reduces to
 *  \f[ Prob(j/i,\lambda)=\sum_{i=1}^4\frac{Z_i Z_i^{-1}}{1 - \psi_i \lambda} \f]
 *  where \f$Z_i \f$ is the matrix of eigenvectors and \f$\psi_i \f$ are the
 *  eigenvalues for the HKY model.
 */
void
update_Q_matrix (phylogeny p, double lambda)
{
	int i,j,k;
	for (i=0; i < 4; i++) {
		for (j=0; j < 4;j++) {
			p->Q[j][i] = 0.;
			for (k=0; k < 4; k++) p->Q[j][i] += (p->z1[k][i] * p->z2[k][j])/(1. + (p->psi[k] * lambda));

		}
	}
}

/*! \brief ln(likelihood) at a node
 *
 *  Calculates \f$ l=log(L)\f$ at a given branch, since every non-root node
 *  has a branch associated with it. \f$ L=\prod_{i=1}^{nsites}L^i\f$ where 
 *  \f[ 
 *  L^i=\sum_{s_1}\sum_{s_2}\pi(s_1)L^i_{up}(s_1)P(s_2/s_1,t)L^i_{down}(s2)
 *  \f]
 *  \return (double) \f$l=\sum_ilog(L^i)\f$
 */

void
ln_likelihood (phylogeny *p, topology t)
{
  int i;
	t->likelihood_proposal = 0.;

  update_likelihood (p, t->root->left, t);

	for (i=t->first_segment; i <= t->last_segment; i++) {
    up   = p[i]->l[0]->u_current->lk;
    down = p[i]->l[t->root->left->id]->d_current->next->lk;

		p[i]->likelihood_proposal = 0.;
		for (site = 0; site < p[i]->nchar; site++) {
			LikSite = 0.;
			for (s1 = 0; s1 < 4; s1++) for (s2 = 0; s2 < 4; s2++)
				LikSite+= p[i]->pi[s1] * up[site][s1] * p[i]->Q[s1][s2] * down[site][s2];
			p[i]->likelihood_proposal	+= ( p[i]->site_weight[site] * log (LikSite) );
		}
		t->likelihood_proposal += p[i]->likelihood_proposal;
	}
}

void
update_likelihood (phylogeny *p, topol_node this, topology t)
{
	int i;
	if (this->left->internal)  update_likelihood_moved_branches (p, this->left, t);
	if (this->right->internal) update_likelihood_moved_branches (p, this->right, t);

	for (i=t->first_segment; i <= t->last_segment; i++) {
    left = p[i]->l[this->left->id]->d_current->next->lk;
    right = p[i]->l[this->right->id]->d_current->next->lk;

    for (site = 0; site < p[i]->nchar; site++) for (s1 = 0; s1 < 4; s1++) {
      lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p[i]->Q[s1][s2] * left[site][s2];
				lkr += p[i]->Q[s1][s2] * right[site][s2];
			} //for (s2 < 4)
			p[i]->l[this->id]->d_current->next->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)
  }
}

void 
accept_likelihood (phylogeny *p, topology t)
{
	int i;

  for (i=t->first_segment; i <= t->last_segment; i++) {
    p[i]->likelihood_current = p[i]->likelihood_proposal;
  }
  t->likelihood_current = t->likelihood_proposal;

  mark_d_done (p, t->root->left, t);

	if (t->root->left->left->internal) 
		update_accepted_likelihood (p, t->root->left->left, t->root, t);
	if (t->root->left->right->internal) 
		update_accepted_likelihood (p, t->root->left->right, t->root, t);
}	

void
update_accepted_likelihood (phylogeny *p, topol_node this, topol_node upnode, topology t)
{
	int i;

	for (i=t->first_segment; i <= t->last_segment; i++) {
		
		up     = p[i]->l[upnode->id]->u_current->lk;
		sister = p[i]->l[this->sister->id]->d_current->lk;

		for (site = 0; site < p[i]->nchar; site++) for (s1 = 0; s1 < 4; s1++) {
			lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p[i]->Q[s1][s2] * up[site][s2];
				lkr += p[i]->Q[s1][s2] * sister[site][s2];
			} //for (s2 < 4)
			p[i]->l[this->id]->u_current->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)

	}

  if (this->left->internal)  update_accepted_likelihood (p, this->left, this, t);
  if (this->right->internal) update_accepted_likelihood (p, this->right, this, t);
}

void
ln_likelihood_moved_branches (phylogeny *p, topology t)
{
	int i;
	topol_node ingroup = t->root->left;
	t->likelihood_proposal = 0.;


	while ((ingroup->internal) && (ingroup->u_done)) {
		if (!ingroup->left->d_done) ingroup = ingroup->left;
		else  ingroup = ingroup->right;
	}
	if (ingroup != t->root->left) ingroup = ingroup->up;


  if (!(ingroup->d_done))
    update_likelihood_moved_branches (p, ingroup, t);
	else {
		for (i=t->first_segment; i <= t->last_segment; i++) {
			p[i]->likelihood_proposal = p[i]->likelihood_current;
			t->likelihood_proposal += p[i]->likelihood_proposal;
		}
		return;
	}


	for (i=t->first_segment; i <= t->last_segment; i++) {
		if (ingroup == t->root->left) 
			up = p[i]->l[0]->u_current->lk;
		else 
			up = p[i]->l[ingroup->id]->u_current->lk;

    if (ingroup->d_done)
      down = p[i]->l[ingroup->id]->d_current->lk;
    else
      down = p[i]->l[ingroup->id]->d_current->next->lk;

		p[i]->likelihood_proposal = 0.;
		for (site = 0; site < p[i]->nchar; site++) {
			LikSite = 0.;
			for (s1 = 0; s1 < 4; s1++) for (s2 = 0; s2 < 4; s2++)
				LikSite+= p[i]->pi[s1] * up[site][s1] * p[i]->Q[s1][s2] * down[site][s2];
			p[i]->likelihood_proposal	+= ( p[i]->site_weight[site] * log (LikSite) );
		}
		t->likelihood_proposal += p[i]->likelihood_proposal;
	}
}

void
update_likelihood_moved_branches (phylogeny *p, topol_node this, topology t)
{
	int i;
	if ((this->left->internal) && (!(this->left->d_done))) {
		update_likelihood_moved_branches (p, this->left, t);
	}
	if ((this->right->internal) && (!(this->right->d_done))) {
		update_likelihood_moved_branches (p, this->right, t);
	}

	for (i=t->first_segment; i <= t->last_segment; i++) {
	
    if (this->left->d_done) 
      left = p[i]->l[this->left->id]->d_current->lk;
    else 
      left = p[i]->l[this->left->id]->d_current->next->lk;
    if (this->right->d_done) 
      right = p[i]->l[this->right->id]->d_current->lk;
    else 
      right = p[i]->l[this->right->id]->d_current->next->lk;


		for (site = 0; site < p[i]->nchar; site++) for (s1 = 0; s1 < 4; s1++) {
			lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p[i]->Q[s1][s2] * left[site][s2];
				lkr += p[i]->Q[s1][s2] * right[site][s2];
			} //for (s2 < 4)
			p[i]->l[this->id]->d_current->next->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)

	}
}

void 
accept_likelihood_moved_branches (phylogeny *p, topology t)
{
	int i;
	topol_node ingroup = t->root->left;
//	bool debg = false;
//	int ud[t->nnodes], dd[t->nnodes];

	for (i=t->first_segment; i <= t->last_segment; i++) {
		p[i]->likelihood_current = p[i]->likelihood_proposal;
	}
	t->likelihood_current = t->likelihood_proposal;
		
/*	for (i=t->nleaves; (i < t->nnodes); i++) { 
		ud[i] = (int) (t->nodelist[i]->u_done); 
		dd[i] = (int) (t->nodelist[i]->d_done); } */


	while ((ingroup->internal) && (ingroup->u_done)) {
		if (!ingroup->left->d_done) ingroup = ingroup->left;
		else  ingroup = ingroup->right;
	}
	if (ingroup != t->root->left) ingroup = ingroup->up;


	if (!(ingroup->d_done)) {
		if (ingroup->internal) mark_d_done_moved_branches (p, ingroup, t);
		if (ingroup != t->root->left) {
			update_likelihood_moved_branches (p, t->root->left, t);
      mark_d_done_moved_branches (p, t->root->left, t);
		}
	}

	t->root->left->u_done = true;
	if (t->root->left->left->internal) 
		update_accepted_likelihood_moved_branches (p, t->root->left->left, t->root, t);
	if (t->root->left->right->internal) 
		update_accepted_likelihood_moved_branches (p, t->root->left->right, t->root, t);
/*
	for (i=t->nleaves; (!debg) && (i < t->nnodes); i++) if (!t->nodelist[i]->d_done) { debg = true; }
	if (debg) {
		for (i=t->nleaves; (i < t->nnodes); i++) {
			printf ("{%d%d|%d%d| %d %d} ", ud[i], (int)(t->nodelist[i]->u_done), 
							dd[i], (int)(t->nodelist[i]->d_done), t->nodelist[i]->id, t->nodelist[i]->up->id);
		}
		printf ("\n[[%d %d %d]] ", t->root->left->left->id, t->root->left->id, t->root->left->right->id);
	}
 */
}	

void
update_accepted_likelihood_moved_branches (phylogeny *p, topol_node this, topol_node upnode, topology t)
{
	int i;

	for (i=t->first_segment; i <= t->last_segment; i++) {
    p[i]->l[this->id]->u_current = p[i]->l[this->id]->u_current->next;

		up     = p[i]->l[upnode->id]->u_current->lk;
		sister = p[i]->l[this->sister->id]->d_current->lk;

		for (site = 0; site < p[i]->nchar; site++) for (s1 = 0; s1 < 4; s1++) {
			lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p[i]->Q[s1][s2] * up[site][s2];
				lkr += p[i]->Q[s1][s2] * sister[site][s2];
			} //for (s2 < 4)
			p[i]->l[this->id]->u_current->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)

	}

	this->u_done = true;
	if (this->left->internal)  update_accepted_likelihood_moved_branches (p, this->left, this, t);
	if (this->right->internal) update_accepted_likelihood_moved_branches (p, this->right, this, t);
}

void
ln_likelihood_proposal (phylogeny p, topology t)
{
	p->likelihood_proposal = 0.;

	update_likelihood_proposal (p, t->root->left);

	up   = p->l[0]->u_accepted->lk;
	down = p->l[t->root->left->id]->d_accepted->next->lk;

	for (site = 0; site < p->nchar; site++) {
		LikSite = 0.;
		for (s1 = 0; s1 < 4; s1++) for (s2 = 0; s2 < 4; s2++)
			LikSite+= (p->pi[s1] * up[site][s1] * p->Q[s1][s2] * down[site][s2]);
		p->likelihood_proposal += p->site_weight[site] * log (LikSite);
	}
}

void
update_likelihood_proposal (phylogeny p, topol_node this)
{
	if (this->left->internal)  update_likelihood_proposal (p, this->left);
	if (this->right->internal) update_likelihood_proposal (p, this->right);
	
	left  = p->l[this->left->id]->d_accepted->next->lk;
	right = p->l[this->right->id]->d_accepted->next->lk;
	
	for (site = 0; site < p->nchar; site++) {
		for (s1 = 0; s1 < 4; s1++) {
			lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p->Q[s1][s2] * left[site][s2];
				lkr += p->Q[s1][s2] * right[site][s2];
			} //for (s2 < 4)
			p->l[this->id]->d_accepted->next->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)
	} // for (site)
}

void 
accept_likelihood_proposal (phylogeny p, topology t)
{
	int i;
	p->likelihood_accepted = p->likelihood_current = p->likelihood_proposal;
	for (i=p->ntax; i < p->nnodes; i++)
		p->l[i]->d_accepted = p->l[i]->d_accepted->next;

	if (t->root->left->left->internal) 
		update_accepted_likelihood_proposal (p, t->root->left->left, t->root);
	if (t->root->left->right->internal) 
		update_accepted_likelihood_proposal (p, t->root->left->right, t->root);
}	


void
update_accepted_likelihood_proposal (phylogeny p, topol_node this, topol_node upnode)
{
	up     = p->l[upnode->id]->u_accepted->lk;
	sister = p->l[this->sister->id]->d_accepted->lk;
	
	for (site = 0; site < p->nchar; site++) {
		for (s1 = 0; s1 < 4; s1++) {
			lkl = lkr = 0.0;
			for (s2 = 0; s2 < 4; s2++) {
				lkl += p->Q[s1][s2] * up[site][s2];
				lkr += p->Q[s1][s2] * sister[site][s2];
			} //for (s2 < 4)
			p->l[this->id]->u_accepted->lk[site][s1] = lkr * lkl;
		} // for (s1 < 4)
	} // for (site)
	
	if (this->left->internal)  update_accepted_likelihood_proposal (p, this->left, this);
	if (this->right->internal) update_accepted_likelihood_proposal (p, this->right, this);
}

void
mark_d_done (phylogeny *p, topol_node this, topology t)
{
	int i;
	if (this->left->internal)  mark_d_done (p, this->left, t);
	if (this->right->internal) mark_d_done (p, this->right, t);

	for (i=t->first_segment; i <= t->last_segment; i++) {
    p[i]->l[this->id]->d_current = p[i]->l[this->id]->d_current->next;
	}
}

void
mark_d_done_moved_branches (phylogeny *p, topol_node this, topology t)
{
	int i;
  if ((this->left->internal) && !(this->left->d_done)) {
    mark_d_done_moved_branches (p, this->left, t);
  }
  if ((this->right->internal) && !(this->right->d_done)) {
    mark_d_done_moved_branches (p, this->right, t);
  }

  for (i=t->first_segment; i <= t->last_segment; i++) {
    p[i]->l[this->id]->d_current = p[i]->l[this->id]->d_current->next;
  }

  this->d_done = true;
}
