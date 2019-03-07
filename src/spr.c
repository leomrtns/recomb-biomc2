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

/*! \file spr.c 
 * \brief Subtree Prune-Regraft operation functions. 
 *
 * For the purposes of branch swapping, "edge" and "node" are
 * synonyms, if we remember that to every non-root node there is an edge
 * connecting it to its parent.
 */

#include "spr.h"

/* Local variables */

/*! \brief Vector with IDs of elegible or non-elegible nodes for swapping. */
int node_urn[8];

/* Local function prototypes */

/* \brief Randomly choose a node different from root. */
topol_node spr_draw_prune_node (topology T);

/*! \brief Randomly choose an internal node. */
topol_node spr_draw_prune_internal_node (topology T);

/*! \brief Randomly choose a node which is not neighbour to prune node */
topol_node spr_draw_regraft_node_spr (topology T, topol_node prune, int *node_index);

/*! \brief Randomly choose a node which is not neighbour to prune node but share
 * a common neighbour. */
topol_node spr_draw_regraft_node_nni (topology T, topol_node prune);
/*! \brief Flag nodes descending from this topol_node_struct as upper part undone, in pre-order. */
void undo_udone (topol_node this);

bool
spr_apply_multiple_moves_eulertour (topology T, int n_moves)
{
	int inodes = T->nleaves-3, i, j, j2, k, k2, n_tries = 0;
	int node_index[T->nnodes-1], sample = 32;
	long long int euler_tour[n_moves+1][sample][inodes];
	split_partition *spt[inodes];
	topol_node prune, regraft;
	bool is_new = true;

	if (n_moves == 1) {
		prune = spr_draw_prune_node (T);
		regraft = spr_draw_regraft_node_spr (T, prune, node_index);
		spr_apply_move_at_nodes (T, prune, regraft, true);
		return true;
	}
	else { /* we do not check, but n_moves should be > 1 */
		create_topology_bipartition_l (T);
		order_topology_by_bipartition (T->root->left);
		i = 0; postorder_mapping_bipartition (T->root->left, (spt), &i);
		/* original topology does not need sampling */
		for (i=0; i < inodes; i++) euler_tour[0][0][i] = spt[i]->l;

		/* first SPR */
		for (k2=0; k2 < sample; k2++) { 
			prune = spr_draw_prune_node (T);
			regraft = spr_draw_regraft_node_spr (T, prune, node_index);
			spr_apply_move_at_nodes (T, prune, regraft, true);

			create_topology_bipartition_l (T);
			order_topology_by_bipartition (T->root->left);
			i = 0; postorder_mapping_bipartition (T->root->left, (spt), &i);
			for (i=0; i < inodes; i++) euler_tour[1][k2][i] = spt[i]->l;
			/* sample #sample topologies, use the last one */
			if (k2 < (sample-1)) spr_apply_move_at_nodes (T, T->undo_prune, T->undo_regraft, true);
		}

		/* for more than one SPR we need to check for cycles */
		for (j=2; j<=n_moves; j++) {
			n_tries = 0;
			for (k2=0; k2 < sample; k2++) { 
				prune = spr_draw_prune_node (T);
				regraft = spr_draw_regraft_node_spr (T, prune, node_index);
				spr_apply_move_at_nodes (T, prune, regraft, true);

				create_topology_bipartition_l (T);
				order_topology_by_bipartition (T->root->left);
				i = 0; postorder_mapping_bipartition (T->root->left, (spt), &i);
				for (i=0; i < inodes; i++) euler_tour[j][k2][i] = spt[i]->l;
				if (k2 < (sample-1)) spr_apply_move_at_nodes (T, T->undo_prune, T->undo_regraft, true);

				/* euler_tour[i] is always != euler_tour [i-1] */
				for (is_new = false, k=0; (k < inodes) && (!is_new); k++)
					if (euler_tour[j][k2][k] != euler_tour[0][0][k]) is_new = true;
				for (i = 1; (i<j-1) && (is_new); i++)
					for (j2=0; (j2 < sample) && (is_new); j2++) 
						for (is_new = false, k=0; (k < inodes) && (!is_new); k++)
							if (euler_tour[j][k2][k] != euler_tour[i][j2][k]) is_new = true;

				if (!is_new) {
					if (n_tries < 16) k2--; 
					n_tries++;
					/* try again the same distance */
					if (k2 == (sample - 1)) {
						if (n_tries > 16) return false;
						spr_apply_move_at_nodes (T, T->undo_prune, T->undo_regraft, true);
						j--;
					}
				}
			}
		
		} // for (j <= n_moves)
		
		return true;
	} // else 
}

bool
spr_apply_multiple_moves (topology T, int n_moves)
{
	int i, j, n_tries = 0;
	int node_index[T->nnodes-1];
	topol_node prune, regraft;
	topol_node p_list[n_moves], r_list[n_moves];
	bool is_new = true;

	if (n_moves > T->nleaves - 3) n_moves = T->nleaves-3;

	if (n_moves == 1) {
		prune = spr_draw_prune_node (T);
		regraft = spr_draw_regraft_node_spr (T, prune, node_index);
		spr_apply_move_at_nodes (T, prune, regraft, true);
//    spr_apply_move_at_nodes (p, p->undo_prune, p->undo_regraft, 0);
		return true;
	}
	else { /* we do not check, but n_moves should be > 2 */
		/* first SPR */
		prune = spr_draw_prune_node (T);
		regraft = spr_draw_regraft_node_spr (T, prune, node_index);
		spr_apply_move_at_nodes (T, prune, regraft, true);
		p_list[0] = T->undo_prune;
		r_list[0] = T->undo_regraft;

		/* for more than one SPR we need to check for cycles */
		for (j=1; j < n_moves; j++) {
			prune = spr_draw_prune_node (T);
			regraft = spr_draw_regraft_node_spr (T, prune, node_index);

			is_new = true;
			for (i=0; i < j && is_new; i++) {
				if (prune   == p_list[i]) is_new = false;
				if (regraft == r_list[i]) is_new = false;
			}
			if (is_new) {
				spr_apply_move_at_nodes (T, prune, regraft, true);
				p_list[j] = T->undo_prune;
				r_list[j] = T->undo_regraft;
				n_tries = 0;
			}
			else {
				if (n_tries > 32) return false;
				n_tries++;
				j--;
			}
		} // for (j <= n_moves)

		return true;
	} // else 
}

void
spr_apply_random_move_internal (topology T, bool done)
{
	int node_index[T->nnodes-1];
  topol_node prune, regraft;
	prune = spr_draw_prune_internal_node (T);
	regraft = spr_draw_regraft_node_spr (T, prune, node_index);
	spr_apply_move_at_nodes (T, prune, regraft, done);
}

void
spr_apply_random_move_spr (topology T, bool done)
{
	int node_index[T->nnodes-1];
  topol_node prune, regraft;
  prune = spr_draw_prune_node (T);
	regraft = spr_draw_regraft_node_spr (T, prune, node_index);
  spr_apply_move_at_nodes (T, prune, regraft, done);
}

void
spr_apply_random_move_nni (topology T, bool done)
{
  topol_node prune, regraft;
  prune = spr_draw_prune_node (T);
	regraft = spr_draw_regraft_node_nni (T, prune);
  spr_apply_move_at_nodes (T, prune, regraft, done);
}

topol_node
spr_draw_prune_node (topology T)
{
	if (T->nleaves == 4) 
   {
    /* situation where node cannot be pruned (no eligible neighbours). Happens
     * only when there are four seqs */
		int i = biomc2_rng_uniform_int (random_number, 4);
		if ( i == 0) return T->root->left;
		else return T->nodelist[i];
	 }// if (workaround for 4 seqs) 
	
	else {
		int i = biomc2_rng_uniform_int (random_number, T->nnodes - 1);
		return T->nodelist[i+1]; /* nodelist[0] is the root node */
	}
}

topol_node
spr_draw_prune_internal_node (topology T)
{
	if (T->nleaves == 4) return spr_draw_prune_node (T);
	else return T->nodelist[T->nleaves + biomc2_rng_uniform_int (random_number, T->nnodes - T->nleaves)];
}

topol_node
spr_draw_regraft_node_spr (topology T, topol_node prune, int *node_index)
{
	int i;
	int n_regraft, current, ne = 0; /* vector with non-eligible nodes */
	
	node_urn[ne++] = prune->id; /* self */

	/* ingroup root, do not have sister, and root was already removed */
	if (prune->up == T->root) {
		node_urn[ne++] = prune->left->id;
		node_urn[ne++] = prune->right->id;
	}
	/* another node, with up to four neighbours */
	else {
		node_urn[ne++] = prune->up->id;
		if (prune->left)  { node_urn[ne++] = prune->left->id; }
		if (prune->right) { node_urn[ne++] = prune->right->id; }
    /* sister */
    if (prune->up->right == prune) { node_urn[ne++] = prune->up->left->id; }
    else { node_urn[ne++] = prune->up->right->id; }
   }
  qsort (node_urn, ne, sizeof (int), compare_int);

	n_regraft = current = 0;
	/* excludes i = 0 since this is root node */
	for (i = 1; i < T->nnodes; i++) {
		if ((current < ne) && (i == node_urn[current])) current++;
		else node_index[n_regraft++] = i;
	}
  
  i = biomc2_rng_uniform_int (random_number, n_regraft);
  
  return T->nodelist[node_index[i]];
}

topol_node
spr_draw_regraft_node_nni (topology T, topol_node prune)
{
	int i;
	int ne = 0; /* vector with eligible nodes */
	topol_node sister;

	if ((prune->up) && (prune->up->up)) {

		if (prune->up->up->up) {
			node_urn[ne++] = prune->up->up->id;
			if (prune->up->up->left == prune->up) node_urn[ne++] = prune->up->up->right->id;
			else node_urn[ne++] = prune->up->up->left->id;
		}

		if (prune->up->left == prune) sister = prune->up->right;
		else sister = prune->up->left;

		if (sister->left)  node_urn[ne++] = sister->left->id;
		if (sister->right) node_urn[ne++] = sister->right->id;
	}

	if (prune->left) {
		if (prune->left->left)  node_urn[ne++] = prune->left->left->id;
		if (prune->left->right) node_urn[ne++] = prune->left->right->id;
	}
	if (prune->right) {
		if (prune->right->left)  node_urn[ne++] = prune->right->left->id;
		if (prune->right->right) node_urn[ne++] = prune->right->right->id;
	}

	i = biomc2_rng_uniform_int (random_number, ne);

	return T->nodelist[node_urn[i]];
}

void
spr_apply_move_at_nodes (topology tree, topol_node prune, topol_node regraft, bool not_phylo)
{
  if (node1_is_child_of_node2 (regraft, prune)) {
    /*! 
     * The actual SPR move on a topology rooted on a leaf (unrooted topology,
     * for phylogenetic purposes) needs to handle two cases: <b>prune node is in the
     * path from regraft node to the root</b> (prune node is least common ancestor between prune and regraft)
     * and <b>prune node is not in the path from regraft node to root</b> (prune
     * and regraft nodes share a common ancestor). 
     *
     * <b> prune is lca</b>:
     * Algorithm equivalent to rerooting, regraft node climbs up till it finds prune.
     */
    /*! \verbatim
     *
     *                        /prune.up        prune.up\
     * regraft_________ prune/            ==>           \prune___________C
     *          |   |        \                          /        |   |  
     *          A   B         \C                regraft/         A   B
     *
     * \endverbatim
     */
		topol_node r           = regraft, 
							 rup         = regraft->up, 
							 tmp         = rup->up, 
							 newchild    = regraft->up, 
							 prunesister = prune->sister;

    r->up = prune;
    rup->up = prune;
		if (rup->left == r) {
			rup->left = tmp;
			rup->right->sister = tmp; tmp->sister = rup->right;
		}
		else {
			rup->right = tmp;
			rup->left->sister = tmp; tmp->sister = rup->left;
		}
		rup->sister = r; r->sister = rup;

    r = rup;
    rup = tmp;
    tmp = tmp->up;
		while (rup != prune) {
      if (rup->left == r) { /* this won't work for last iteration, see fix [1,2] below */ 
				rup->left = tmp; 
				rup->right->sister = tmp; tmp->sister = rup->right;
			}
			else {
				rup->right = tmp;
				rup->left->sister = tmp; tmp->sister = rup->left;
			}

      rup->up = r;
			r = rup;
			rup = tmp;
			tmp = tmp->up;
		}
		/* rup == prune */
		if (prune->left == r) tmp = prune->right; /* sister */
    else tmp = prune->left;
    tmp->up = r;
		if (r->left == prune) { /* [1] fix commented above */
			r->left = tmp;
			r->right->sister = tmp; tmp->sister = r->right;
		}
		else {
			r->right = tmp;
			r->left->sister = tmp; tmp->sister = r->left;
		}
		prune->sister = prunesister; /* [2] recover original sister */

		/* we can arbitrarily set prune->left = regraft */
		prune->left = regraft;
		prune->right = newchild;

    /* how to undo this move */
    tree->undo_prune = prune;
    tree->undo_regraft = tmp;
  
		if (!not_phylo) {
			if (prune->left->internal)  undo_udone (prune->left); 
			if (prune->right->internal) undo_udone (prune->right);
			if (!(tmp->internal)) tmp = tmp->up;
			while (tmp->up) { tmp->d_done = not_phylo; tmp = tmp->up; }
		}

	}

  else {
    /*! <b>prune is not lca</b>:
     *  Detach the prune subtree and reinsert it just above the regraft node.
     */
    /*! \verbatim
     *  Prune: 
     *
     *  p.left\              /p.up.up                       p.left\                |p.up.up                  
     *         \prune___p.up/                       ==>            \p_______prune  |                         
     *         /            \                                      /               |                         
     * p.right/              \p.up.left || p.up.right      p.right/                |p.up.left || p.up.right  
     
     \endverbatim
     */
    topol_node psister, p = prune;

		if (!not_phylo) {
			prune = prune->up;
			while (prune->up) { prune->d_done = not_phylo; prune = prune->up; }
			prune = regraft->up;
			while ((prune->up) && (prune->d_done)) { prune->d_done = not_phylo; prune = prune->up; }
			/* prune == LCA between original prune and regraft nodes */
			if (prune->left->internal) undo_udone (prune->left);
			if ((prune->right) && (prune->right->internal)) undo_udone (prune->right);
			prune = p;
		}

		psister = prune->sister;
    prune = prune->up;
		psister->sister = prune->sister;
		if (prune->sister) prune->sister->sister = psister;
    psister->up = prune->up;
    if (prune->up->left == prune)
      prune->up->left = psister;
    else
      prune->up->right = psister;
    
    /*! \verbatim
     *  Regraft: 
     *
     *  p.left\                |r.up        p.left\               /prune.up (=r.up.up)
     *         \p_______prune  |      ==>          \p_______prune/ 
     *         /               |                   /             \ 
     * p.right/                |r          p.right/               \r
     
     \endverbatim
     */
    if (prune->left == psister)
      prune->left = regraft;
    else
      prune->right = regraft;

		p->sister = regraft;
		regraft->sister = p;

    prune->up = regraft->up;
		if (regraft->up->left == regraft) {
			regraft->up->left = prune;
			if (regraft->up->right) regraft->up->right->sister = prune;
			prune->sister = regraft->up->right;
		}
		else {
			regraft->up->right = prune;
			regraft->up->left->sister = prune;
			prune->sister = regraft->up->left;
		}
    regraft->up = prune;

    /* how to undo this move */
    tree->undo_prune = p;
    tree->undo_regraft = psister;

	}
//	update_topology_sisters (tree);
}

void
undo_udone (topol_node this)
{
	this->u_done = false;
	if (this->left->internal) undo_udone (this->left);
	if (this->right->internal) undo_udone (this->right);
}
