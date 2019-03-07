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

/*! \file topology_distance.c
 *  \ingroup topology 
 *  \brief Calculation of distances between topologies.
 */

#include "topology.h"

/* local funtion prototypes */

/*! \brief Comparison between split_partition_struct elements based on largest
 * leaf. 
 *
 * Comparison function used by sort(). Results in crescent order (from
 * bipartitions with smallest leaves to bipartitions with largest leaves), with
 * arbitrary ordering of leaves. */
int compare_split_partition_leaf_l (const void *a, const void *b);

/*! \brief Comparison between split_partition_struct elements based on smallest
 * subtree. 
 *
 * Comparison function used by sort(). Results in crescent order (from
 * bipartitions with fewer leaves to bipartitions with more leaves), with
 * ties broken by largest leaf */
int compare_split_partition_l (const void *a, const void *b);

/*! \brief Stores bipartition information in post-order traversal where heavy
 * child (topol_node_struct with more leaves) is stored before light child. */
void heavychild_mapping_bipartition (topol_node this, split_partition **u, int *n);

/*! \brief Calculates approximate \f$d_{SPR}\f$ based on split_space_struct
 * information
 *
 * This function does the actual \f$d_{SPR}\f$ calculation iteratively untill
 * all bipartitions are in agreement or \f$max\f$ iterations are achieved. 
 * \param[in] split split_space_struct with bipartition vectors from two
 * topoloiges
 * \param[in] max maximum number of iterations (usually \f$N-3\f$ where \f$N\f$
 * is the number of leaves
 * \return \f$d_{SPR}\f$ estimate
 */
int dSPR_partition_l (split_space split, int *max);

/*! \brief Auxiliary function to create_topology_bipartition(). */
void create_subtree_bipartition_l (topol_node this, int *n_leaves);

int
compare_split_partition_leaf_l (const void *a, const void *b)
{
	if ( (*(split_partition **)a)->l > (*(split_partition **)b)->l ) return 1;
	if ( (*(split_partition **)a)->l < (*(split_partition **)b)->l ) return -1;
	return 0;
}

int
compare_split_partition_l (const void *a, const void *b)
{
	int r;
	if (( r = (*(split_partition **)a)->n - (*(split_partition **)b)->n))
		return r;
	if ( (*(split_partition **)a)->l > (*(split_partition **)b)->l ) return 1;
	if ( (*(split_partition **)a)->l < (*(split_partition **)b)->l ) return -1;
	return 0;
}

void
create_topology_bipartition_l (topology tree)
{
  create_subtree_bipartition_l (tree->root->left, &(tree->nleaves));
}

void
create_subtree_bipartition_l (topol_node this, int *n_leaves)
{
	if (this->left->internal) create_subtree_bipartition_l (this->left, n_leaves);
	if (this->right->internal) create_subtree_bipartition_l (this->right, n_leaves);

	split_OR_l (&(this->bipartition), &(this->left->bipartition), &(this->right->bipartition));
}

void
order_topology_by_bipartition (topol_node this)
{
	if (this->left->internal) order_topology_by_bipartition (this->left); 
	if (this->right->internal) order_topology_by_bipartition (this->right); 

	if (this->left->bipartition.l > this->right->bipartition.l) {
		topol_node tmp = this->left;
		this->left = this->right;
		this->right = tmp;
	}
}

void
postorder_mapping_bipartition (topol_node this, split_partition **u, int *n)
{
	if (this->left->internal) postorder_mapping_bipartition (this->left, u, n); 
	if (this->right->internal) postorder_mapping_bipartition (this->right, u, n);
	if (!this->up->up) return;

	u[(*n)++] = &(this->bipartition);
}

void
heavychild_mapping_bipartition (topol_node this, split_partition **u, int *n)
{
  if (this->left->bipartition.n > this->right->bipartition.n) {
    if (this->left->internal) heavychild_mapping_bipartition (this->left, u, n); 
    if (this->right->internal) heavychild_mapping_bipartition (this->right, u, n);
  }
  else if (this->left->bipartition.n < this->right->bipartition.n) {
    if (this->right->internal) heavychild_mapping_bipartition (this->right, u, n);
    if (this->left->internal) heavychild_mapping_bipartition (this->left, u, n); 
  }
  /* balanced subtree ( size(left) == size(right) ) */
  else if (this->left->bipartition.l > this->right->bipartition.l) {
    if (this->left->internal) heavychild_mapping_bipartition (this->left, u, n); 
    if (this->right->internal) heavychild_mapping_bipartition (this->right, u, n);
  }
  else { 
    if (this->right->internal) heavychild_mapping_bipartition (this->right, u, n);
    if (this->left->internal) heavychild_mapping_bipartition (this->left, u, n); 
  }

	if (!this->up->up) return;

	u[(*n)++] = &(this->bipartition);
}

void
euler_semitour_mapping (topol_node this, int *u, int *n)
{
  u[(*n)++] = this->id;
	if (this->left)  euler_semitour_mapping (this->left, u, n); 
	if (this->right) euler_semitour_mapping (this->right, u, n);
}

void
update_topology_from_topology_split (topology t1, topology t2, split_space split)
{
	int i;

	/* t1 = from_tree, t2 = to_tree and u_node, d_done annotation is done on t1 */

  initialize_split_space_l (split);

	create_topology_bipartition_l (t1);
	create_topology_bipartition_l (t2);
  
	i = 0; postorder_mapping_bipartition (t1->root->left, split->u1, &i);
	i = 0; postorder_mapping_bipartition (t2->root->left, split->u2, &i);
	
	/* minimize distance from leaves */
	split_minimize_initial_subtree_size_l (split->u1, &(split->n_u1), split);
	split_minimize_initial_subtree_size_l (split->u2, &(split->n_u2), split);
	/* remove identical edges and map IDs between t1 and t2 */
	split_create_agreement_list_mapping (split);
	/* ingroup root node is a trivial bipartition */
	t1->root->left->map_id = t2->root->left->id;
	t2->root->left->map_id = t1->root->left->id;
	t1->root->left->u_done = true; 

	if (split->n_agree == split->inodes) { /* same topology */
		for (i=t1->nleaves; i < t1->nnodes; i++)
			t1->nodelist[i]->u_done = t1->nodelist[i]->d_done = true;
	}
	else {
		/* subtree collapsing with annotation and IDs from topology t1 (from_tree) */ 
		split_compress_subtree_agreement_topol (split);
	}
}

bool
dSPR_is_zero_l (topology t1, topology t2, split_space split)
{
	int i;

  initialize_split_space_l (split);

	create_topology_bipartition_l (t1);
	create_topology_bipartition_l (t2);
  
	i = 0; heavychild_mapping_bipartition (t1->root->left, split->u1, &i);
	i = 0; heavychild_mapping_bipartition (t2->root->left, split->u2, &i);
	
  // DEBUG:
  //for (i=0; i < split->n_u1; i++) split_print_binary_bipartition_l ((split->u1[i]), split); printf (" [t1] \n");
  //for (i=0; i < split->n_u2; i++) split_print_binary_bipartition_l ((split->u2[i]), split); printf (" [t2] \n");
  
  for (i=0; i < split->n_u1; i++) if (split->u1[i]->l != split->u2[i]->l) return false;

  return true; /* else */
}

int
dSPR_topology_l (topology t1, topology t2, split_space split, int *max)
{
	bool isZero = dSPR_is_zero_l (t1, t2, split); /* initialize vectors */

  if (isZero) return 0;  /* irrespective of (*max) the answer is zero */
  if (!(*max)) return 1; /* (*max) == zero is equiv to !isZero, and we know it's not zero */
  return dSPR_partition_l (split, max); /* else (the topols are diff and we should quantify it) */
}

int 
dSPR_partition_l (split_space split, int *max)
{
	int i; 
	bool mismatch_edges = true;
  bool firstpass = true, verbose = 0; /* DEBUG option so that we can see the steps of d_SPR calc */

	if (verbose > 3) {
		for (i=0; i < split->n_u1; i++) split_print_binary_bipartition_l ((split->u1[i]), split); printf (" [t1] \n");
		for (i=0; i < split->n_u2; i++) split_print_binary_bipartition_l ((split->u2[i]), split); printf (" [t2] \n\n");
	}

	/* minimize distance from leaves */
	split_minimize_initial_subtree_size_l (split->u1, &(split->n_u1), split);
	split_minimize_initial_subtree_size_l (split->u2, &(split->n_u2), split);

  split->rfdistance = (int) (((double)(split->n_u1 + split->n_u2))/2.);
	mismatch_edges = ((split->n_u1 > 0) && (split->n_u2 > 0));

  while (mismatch_edges && (split->dSPR <= (*max))) {
		/* remove identical edges*/
		split_create_agreement_list_l (split);

		/* remove agreement edges, for when the identical edge appears more than once (e.g. in minimized tree, after a few rounds)*/
		split_remove_agreement_edges_l (split->u1, &(split->n_u1), split); 
		split_remove_agreement_edges_l (split->u2, &(split->n_u2), split); 

		/* subtrees with two leaves are replaced by a new "leaf" (subtree collapse) */
		split_compress_subtree_agreement_l (split);
		
		/* leaf compression will change the sizes */
//		split_minimize_subtree_size_l (split->u1, &(split->n_u1), split);
//		split_minimize_subtree_size_l (split->u2, &(split->n_u2), split);
		/* sort trees (speed up pairwise comparisons) */ 
//		qsort ((split->u1), split->n_u1, sizeof (split_partition*), compare_split_partition_l);
//		qsort ((split->u2), split->n_u2, sizeof (split_partition*), compare_split_partition_l);

		if (verbose > 1) {
			if (verbose > 2) {
				for (i=0; i < split->n_u1; i++) split_print_binary_bipartition_l ((split->u1[i]), split); printf (" [t1] \n");
				for (i=0; i < split->n_u2; i++) split_print_binary_bipartition_l ((split->u2[i]), split); printf (" [t2] \n\n");
			}
			for (i=0; i < split->n_agree; i++) split_print_binary_bipartition_l (split->agree[i], split); printf (" [agree] \n");
			for (i=0; i < split->n_leaves; i++) printf ("%*d.", 2, i+1);
			printf ("\n");
			for (i=0; i < split->n_leaves; i++) printf ("%*d.", 2,split->leaf_idx[i]);
			printf ("[original leaves]\n\n");
		}

    if (firstpass) { split->rfdistance = (int) (((double)(split->n_u1 + split->n_u2))/2.); firstpass = false; }
		mismatch_edges = ((split->n_u1 > 0) && (split->n_u2 > 0));
		if (!mismatch_edges)	return split->dSPR;

		/* compare possible disagreements */
		split_create_disagreement_list_l (split);

		if (split->n_disagree > 1) {
			qsort ((split->d_ptr), split->n_disagree, sizeof (split_partition*), compare_split_partition_l);
			split_remove_duplicate_disagreement_l (split);
		}

		if (split->d_ptr[0]->n > 1) split_find_smallest_disagreement_subtree_l (split);
		else  split->best_prune = split->d_ptr[0]; 

		if (verbose > 1) {
			if (verbose > 2) for (i=0; i < split->n_disagree; i++) { 
				split_print_binary_bipartition_l ((split->d_ptr[i]), split); printf ("\n");
			}
			printf ("\t{ %d }\t best: ", split->n_disagree); 
			split_print_binary_bipartition_l ((split->best_prune), split); printf ("\n");
		}

		split_remove_pruned_subtree_l (split->best_prune, split);
//		split->dSPR += split->best_prune->n; 
		split->dSPR++;
		
		if (verbose > 1) {
			for (i=0; i < split->n_leaves; i++) printf ("%*d.", 2, split->leaf_recomb[i]);
			printf ("[recomb leaves]\n\n");
		}

		/* minimize distance from leaves and remove single leaves */
		split_minimize_subtree_size_l (split->u1, &(split->n_u1), split);
		split_minimize_subtree_size_l (split->u2, &(split->n_u2), split);
		split_minimize_subtree_size_l (split->agree, &(split->n_agree), split);

		mismatch_edges = ((split->n_u1 > 0) && (split->n_u2 > 0));

	} // while (mismatch_edges && SPR <= max)
	return split->dSPR;
}

