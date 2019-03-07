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

/*! \file topology_split_partition.c
 *  \ingroup topology 
 *  \brief split_partition_struct (edge bipartition) calculation.
 */

#include "topology.h"

/*! \brief Replace information of redundant leaf by information of largest leaf
 * for split_space_struct. */
void split_space_replace_bit (int *to, int *from, split_space split);

/*! \brief Replace information of redundant leaf by information of largest leaf
 * for one element of split_space_struct. */
void split_bipartition_replace_bit (split_partition *u, int *from, long long int *masktrue, long long int *maskfalse);

/*! \brief Maps original leaves to compressed leaves into split_space::leaf_idx . */
void store_compression_info_real (int *to, int *from, split_space split);

/*! \brief Does nothing (if \f$ d_{SPR} \f$ extension is not used) */
void store_compression_info_dummy (int *to, int *from, split_space split);

/*! \brief Uses split_space::best_prune and split_space::leaf_idx to store
 * pruned subtree (original leaves) information into split_space::leaf_recomb . */
void extend_subtree_info_real (int *index, int *idx_n, split_space split);

/*! \brief Does nothing (if \f$ d_{SPR} \f$ extension is not used) */
void extend_subtree_info_dummy (int *index, int *idx_n, split_space split);

split_space
new_split_space (int n_leaves)
{
	split_space split;

	split = 
	(split_space) biomc2_malloc (sizeof (struct split_space_struct));

	split->n_leaves = n_leaves;
	split->inodes = n_leaves - 3;
	split->n_active = n_leaves;
  split->dSPR = split->rfdistance = 0;

  /* pointers to topologies bipartitions */
	split->n_u1 = split->n_u2 = split->inodes;
	split->u1 = 
	(split_partition**) biomc2_malloc (split->n_u1 * sizeof (split_partition*));
	split->u2 = 
	(split_partition**) biomc2_malloc (split->n_u2 * sizeof (split_partition*));

  /* index of original leaves mapped to compressed leaves */
	split->leaf_idx =
	(int*) biomc2_malloc (split->n_active * sizeof (int));
	/* leaves involved in recombination */
	split->leaf_recomb =
	(int*) biomc2_malloc (split->n_active * sizeof (int));

  /* list of agreement edges (equal splits on both topologies) */
	split->n_agree = 0;
	split->agree =
  (split_partition**) biomc2_malloc (split->inodes * sizeof (split_partition*));

  /* list of disagreement edges */
  split->n_disagree = 0;
  split->dis1 =
  (split_partition*) biomc2_malloc (split->inodes * split->inodes * sizeof (split_partition));
  split->dis2 =
  (split_partition*) biomc2_malloc (split->inodes * split->inodes * sizeof (split_partition));
  split->d_ptr =
  (split_partition**) biomc2_malloc (split->inodes * split->inodes * sizeof (split_partition*));

	/* default is to neither map original leaves nor extend algo to give pruned subtrees */
	setup_split_dont_store_recomb_info (split);
	return split;
}

void
initialize_split_space_l (split_space split)
{
	int i;

	split->n_active = split->n_leaves;
	split->n_agree = split->n_disagree = 0;
	split->n_u1 = split->n_u2 = split->inodes;
	split->dSPR = split->rfdistance = 0;
  split->active_l = 0LL;

  for (i=0; i < split->n_active; i++) { 
		split->active_l |= (1LL << i); 
	}
}

void
setup_split_store_recomb_info (split_space split)
{
	int i;
	for (i=0; i < split->n_leaves; i++) { 
		split->leaf_idx[i] = i;
		split->leaf_recomb[i] = 0;
	}
	split->extend_subtree_info = &extend_subtree_info_real;
	split->store_compression_info = &store_compression_info_real;
}

void
setup_split_dont_store_recomb_info (split_space split)
{
	split->extend_subtree_info = &extend_subtree_info_dummy;
	split->store_compression_info = &store_compression_info_dummy;
}

void
del_split_space (split_space split)
{
	if (split) {
		if (split->dis2)  free (split->dis2);
		if (split->dis1)  free (split->dis1);
		if (split->u2)    free (split->u2);
		if (split->u1)    free (split->u1);
		if (split->d_ptr) free (split->d_ptr);
		if (split->agree) free (split->agree);
		if (split->leaf_idx)    free (split->leaf_idx);
		if (split->leaf_recomb) free (split->leaf_recomb);
		free (split);
	}
}

void
split_initialize_l (split_partition *bpt, int label_id)
{
  bpt->l = (1LL << label_id);
  bpt->n = 1;
}

void 
split_OR_l (split_partition *res, const split_partition *b1, const split_partition *b2)
{
  res->l = b1->l | b2->l;
  res->n = b1->n + b2->n;
}

void
split_minimize_initial_subtree_size_l (split_partition **u, int *n, split_space split)
{
  int i, k;
  long long int smss_up;

  for (i=0; i < *n; i++) { 
		smss_up = ~u[i]->l;
		k = split->n_active - u[i]->n;

    /* smaller subtree */
    if (k < u[i]->n) { u[i]->n = k; u[i]->l = smss_up; }
		/* break ties by subtree not including root (id = 0)  */
    else if ( (k == u[i]->n) && (u[i]->l & 1LL) ) { u[i]->n = k; u[i]->l = smss_up; }
  }
}

void
split_minimize_subtree_size_l (split_partition **u, int *n, split_space split)
{
  int i, k;
  long long int smss_up;

  for (i=0; i < *n; i++) { 
    smss_up = split->active_l & ~u[i]->l;
    u[i]->l = split->active_l &  u[i]->l;

    split_bipartition_size_l (&(u[i]->n), &(u[i]->l), split);
		if (u[i]->n > 1) k = split->n_active - u[i]->n;
		else k = 0;

    /* smaller subtree */
    if ((k < u[i]->n) && (k > 1)) { u[i]->n = k; u[i]->l = smss_up; }
    /* break ties by subtree that excludes root leaf */
    else if ( (k == u[i]->n) && (k > 1) && (u[i]->l & 1LL) ) { u[i]->n = k; u[i]->l = smss_up; }
    else if ((u[i]->n <= 1) || (k <= 1)) { u[i--] = u[--(*n)]; }
  }
}

void
split_create_agreement_list_l (split_space split)
{
	int i, j;

	for (i=0; i < split->n_u1; i++) 
		for (j=0; (j < split->n_u2); j++) {
			if (split->u1[i]->n == split->u2[j]->n) {
				if (split->u1[i]->l == split->u2[j]->l) {
					split->agree[split->n_agree++] = split->u1[i];
					split->u1[i] = split->u1[--split->n_u1];
					split->u2[j] = split->u2[--split->n_u2];
					i--; j = split->n_u2;
				}
			}
		}
	
}

void
split_create_agreement_list_mapping (split_space split)
{
	int i, j;

	for (i=0; i < split->n_u1; i++) 
		for (j=0; (j < split->n_u2); j++) {
			if (split->u1[i]->n == split->u2[j]->n) {
				if (split->u1[i]->l == split->u2[j]->l) {
					split->u1[i]->node->map_id = split->u2[j]->node->id;
					split->u2[j]->node->map_id = split->u1[i]->node->id;
					split->agree[split->n_agree++] = split->u1[i];
					split->u1[i] = split->u1[--split->n_u1];
					split->u2[j] = split->u2[--split->n_u2];
					i--; j = split->n_u2;
				}
			}
		}
	/* don't have a match but the ID is the same (no swap necessary) */
	for (i=0; i < split->n_u1; i++) 
		for (j=0; (j < split->n_u2); j++) {
			if (split->u1[i]->node->id == split->u2[j]->node->id) {
				split->u1[i]->node->map_id = split->u2[j]->node->id;
				split->u2[j]->node->map_id = split->u1[i]->node->id;
				split->u1[i] = split->u1[--split->n_u1];
				split->u2[j] = split->u2[--split->n_u2];
				i--; j = split->n_u2;
			}
		} 
	if (split->n_u1 != split->n_u2) biomc2_error ( "problem matching IDs between t1 and t2");
	/* don't have a match nor same ID */
	for (i=0; i < split->n_u1; i++) {
		split->u1[i]->node->map_id = split->u2[i]->node->id;
		split->u2[i]->node->map_id = split->u1[i]->node->id;
	}
}

void
split_remove_agreement_edges_l (split_partition **u, int *n, split_space split)
{
  int i, j;

	for (i=0; (i < (*n)); i++) 
		for (j=0; j < split->n_agree; j++)
			if (u[i]->n == split->agree[j]->n) {
				if (u[i]->l == split->agree[j]->l) { u[i--] = u[--(*n)]; j = split->n_agree; }
			}
}

void
split_compress_subtree_agreement_topol (split_space split)
{
	int i, j, k, leafpair[2];

	for (i=0; i < split->n_agree; i++) if (split->agree[i]->n == 2) {
		k = 0;
		for (j=0; j < split->n_active; j++) if ( ((split->agree[i]->l >> j) & 1LL) ) leafpair[k++] = j;
		if (leafpair[0]) { /* root not included */
			if ((split->agree[i]->node->left->d_done) && (split->agree[i]->node->right->d_done)) 
				split->agree[i]->node->d_done = true;
		}
		else { /* root belongs to agreement */
			if ((split->agree[i]->node->up->u_done) && (split->agree[i]->node->sister->d_done)) 
				split->agree[i]->node->u_done = true;
		}

		// collapse info of leaves pair[0] and pair[1] into new leaf at position pair[0]
		split_space_replace_bit (&(leafpair[0]), &(leafpair[1]), split);
		
		// since position pair[1] is now redundant, reduce bipartition space
		split->n_active--; split->active_l = 0LL;
		for (i=0; i < split->n_active; i++) split->active_l |= (1LL << i);
		if (leafpair[1] < split->n_active) 
			split_space_replace_bit (&(leafpair[1]), &(split->n_active), split);

		// recalculate sizes
		split_minimize_subtree_size_l (split->agree, &(split->n_agree), split);
		i = -1; 
	}
}

void
split_compress_subtree_agreement_l (split_space split)
{
	int i, j, k, leafpair[2];

	for (i=0; i < split->n_agree; i++) if (split->agree[i]->n == 2) {
		k = 0;
		for (j=0; j < split->n_active; j++) if ( ((split->agree[i]->l >> j) & 1LL) ) leafpair[k++] = j;

		// collapse info of leaves pair[0] and pair[1] into new leaf at position pair[0]
		split_space_replace_bit (&(leafpair[0]), &(leafpair[1]), split); /* 2013.03: I think this is redundant -- they already have same info */
		
		// since position pair[1] is now redundant, reduce bipartition space
		split->n_active--; split->active_l = 0LL;
		for (i=0; i < split->n_active; i++) split->active_l |= (1LL << i);
		if (leafpair[1] < split->n_active) 
			split_space_replace_bit (&(leafpair[1]), &(split->n_active), split);

		// recalculate sizes
		split_minimize_subtree_size_l (split->agree, &(split->n_agree), split);
		i = -1; 
	}

}

void
split_create_disagreement_list_l (split_space split)
{
	int i, j, k = split->n_u2;

	split->n_disagree = 0;
	for (i=0; (i < split->n_u1) && (k); i++) {
		for (j=0; (j < k); j++) {
			
			split->dis1[split->n_disagree].l = split->active_l & (split->u1[i]->l ^ split->u2[j]->l); 
			split_bipartition_size_l (&(split->dis1[split->n_disagree].n), &(split->dis1[split->n_disagree].l), split);

/*			if (split->dis1[split->n_disagree].n == 1) {
				split->d_ptr[0] = &(split->dis1[split->n_disagree]);
				split->n_disagree = 1;
				return;
			}
*/
			split->dis2[split->n_disagree].l = split->active_l & (split->u1[i]->l ^ ~split->u2[j]->l); 
			split->dis2[split->n_disagree].n = split->n_active - split->dis1[split->n_disagree].n;
/*
			if (split->dis2[split->n_disagree].n == 1) {
				split->d_ptr[0] = &(split->dis2[split->n_disagree]);
				split->n_disagree = 1;
				return;
			}
*/
			if (split->dis1[split->n_disagree].n > split->dis2[split->n_disagree].n) {
				split->d_ptr[split->n_disagree]  = &(split->dis2[split->n_disagree]);
			}
			else {
				split->d_ptr[split->n_disagree]  = &(split->dis1[split->n_disagree]);
			}

			if (split->d_ptr[split->n_disagree]->n <= split->n_active/2) split->n_disagree++;
		}
	}
}

void
split_remove_duplicate_disagreement_l (split_space split)
{
	int i,j;

	for (i = split->n_disagree - 1; i >= 1; i--) 
		if (split->d_ptr[i]->l == split->d_ptr[i-1]->l) {
			for (j = i; j < split->n_disagree-1; j++) split->d_ptr[j] = split->d_ptr[j+1];
			split->n_disagree--;
		}
}

void
split_find_smallest_disagreement_subtree_l (split_space split)
{
	int i, j;
	long long int dis;

	split->best_prune = split->d_ptr[0];
	
	for (i=0; i < split->n_disagree; i++) {
		for (j=0; j < split->n_agree; j++) {
			if ((split->d_ptr[i]->n == split->agree[j]->n) || (split->d_ptr[i]->n == (split->n_active - split->agree[j]->n))) {
				dis = ( split->active_l & (split->d_ptr[i]->l ^ split->agree[j]->l) );
				if (!dis) { split->best_prune = split->d_ptr[i]; return; }
				dis = ( split->active_l & (split->d_ptr[i]->l ^ ~split->agree[j]->l) );
				if (!dis) { 
					split->d_ptr[i]->l = split->active_l & ~split->d_ptr[i]->l;
					split->d_ptr[i]->n = split->n_active - split->d_ptr[i]->n;
					split->best_prune = split->d_ptr[i]; return; }
			}
		}
	}
	/* only if there is no match */
//	if (split->best_prune->n > 2) split_minimize_disagreement_l (split);
}

// UNUSED
void
split_minimize_disagreement_l (split_space split)
{
	int i, size, best_size;
	long long int dis, best_dis;

	best_size = split->best_prune->n; best_dis = split->best_prune->l; 
	for (i=0; (i < split->n_agree) && (best_size > 1) ; i++) if (split->agree[i]->n <= split->best_prune->n) {
		dis = ( split->active_l & (split->agree[i]->l ^ split->best_prune->l) );
		split_bipartition_size_l (&size, &dis, split);
		
		if (size < best_size) { best_size = size; best_dis  = dis; }
	}
	split->best_prune->l = best_dis;
	split->best_prune->n = best_size;
}

void
split_remove_pruned_subtree_l (split_partition *prune, split_space split)
{
	int i, j, index[prune->n], idx_n = 0;

	for (i=0; i < (split->n_active); i++) 
		if ( ((prune->l >> i) & 1LL) ) index[idx_n++] = i; //at the end, idx_n = prune->n 

//	if (idx_n > 1) qsort (index, idx_n, sizeof (int), compare_int);

	// store recomb subtrees (extension to d_SPR) 
	split->extend_subtree_info (index, &idx_n, split);

	for (j = prune->n - 1, i = split->n_active - 1, idx_n = 0; i >= (split->n_active - prune->n); i--) {
		if (index[idx_n] >= (split->n_active - prune->n)) i = -1;
		else {
			if ( i == index[j]) j--;
			else { split_space_replace_bit (&(index[idx_n]), &i, split); idx_n++;}
		}
	}
 
	split->n_active -= prune->n;
	split->active_l = 0LL;
	for (i=0; i < split->n_active; i++) split->active_l |= (1LL << i);
}

void 
split_space_replace_bit (int *to, int *from, split_space split)
{
	int i;
	long long int mtru = ( 1LL << *to), mfal = ~mtru;

	// mapping between original leaves and compressed info
	split->store_compression_info (to, from, split);

	// replace bits 
	for (i=0; i < split->n_u1; i++) split_bipartition_replace_bit (split->u1[i], from, &mtru, &mfal);
	for (i=0; i < split->n_u2; i++) split_bipartition_replace_bit (split->u2[i], from, &mtru, &mfal);
	for (i=0; i < split->n_agree; i++) split_bipartition_replace_bit (split->agree[i], from, &mtru, &mfal);
}

void
split_bipartition_replace_bit (split_partition *u, int *from, long long int *masktrue, long long int *maskfalse)
{
	// masktrue = 0100 , maskfalse = 1011
	if ( ((u->l >> *from) & 1LL) ) u->l |= *masktrue;
	else u->l &= *maskfalse; 
}

void
split_print_binary_bipartition_l (const split_partition *u, split_space split)
{
  int i;
	for (i=0; i < split->n_active; i++)
		printf ("%d", (int)((u->l >> i) & 1LL) );
	printf ("[%d] ", u->n);
}

void
split_bipartition_size_l (int *count, const long long int  *l, split_space split)
{
	int i;
	(*count) = 0;
	for (i=0; (i < split->n_active) ; i++) (*count) += (((*l) >> i) & 1LL);
}

void
extend_subtree_info_real (int *index, int *idx_n, split_space split)
{
	int i, j;
//	for (i=0; i < (*idx_n); i++) printf ("[ %d ]", index[i]);
//	printf ("\n");
	for (i=0; i < (*idx_n); i++) for (j=0; j < split->n_leaves; j++) {
		if (split->leaf_idx[j] == index[i]) split->leaf_recomb[j] = 1L;
	}
}

void
extend_subtree_info_dummy (int *index, int *idx_n, split_space split)
{ (void) index; (void) idx_n; (void) split; return; }


void
store_compression_info_real (int *to, int *from, split_space split)
{ 
	/* before: 0 1 2 3 4 5
	 * after:  0 0 2 3 4 5 (assuming to=0 and from=1) */
	int i;
	for (i=0; i < split->n_leaves; i++) { 
		if (split->leaf_idx[i] == (*from)) split->leaf_idx[i] = (*to);
	}
  /* note from 2013.05.19: leaf_idx[] has the compressed index of each leaf, so that after a compression a 
   * leaf_idx[]= { 0 0 2 0 4 5} means that leaves 0, 1 and 3 are (a common subtree?) represented by position 0 of the bitstring */
}

void
store_compression_info_dummy (int *to, int *from, split_space split)
{ (void) to; (void) from; (void) split;  return; }
