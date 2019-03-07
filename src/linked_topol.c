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

/*! \file linked_topol.c 
 *  \brief 
 */

#include "linked_topol.h"

void shuffle_topology (topology t, phylogeny *segment, linked_topol tp, split_space split);
void insert_topology_right (topology t, linked_topol tp);
void insert_topology_left  (topology t, linked_topol tp);

linked_topol
new_linked_topol (int nsegments, int nleaves)
{
	int i;
	linked_topol tp;

	tp = (linked_topol) biomc2_malloc (sizeof (struct linked_topol_struct));

	/* BKPSIZE are the last topologies in list; used as pivot (store original info) */ 
	tp->size = nsegments + BKPSIZE + 1;
	
	tp->nSPR = tp->nCOP = 0;
	tp->tree = (topology*) biomc2_malloc (tp->size * sizeof (topology));

	for (i=0; i < tp->size; i++) {
		tp->tree[i] = new_topology (nleaves);
		tp->tree[i]->id = i;
	}
	tp->tree[0]->first_segment = 0;
	tp->tree[0]->last_segment  = nsegments-1;
	tp->first = tp->last = tp->tree[0];
	tp->tree[0]->next = tp->tree[0]->prev = NULL;

	tp->tree[1]->next = tp->tree[2];
	for (i=2; i < tp->size-1; i++) {
		tp->tree[i]->prev = tp->tree[i-1];
		tp->tree[i]->next = tp->tree[i+1];
	}
	tp->tree[tp->size-1]->prev = tp->tree[tp->size-2];
	tp->tree[1]->prev = tp->tree[tp->size-1]->next = NULL;

	tp->avail = tp->tree[1];

	for (i=0; i < BKPSIZE; i++) { tp->bkp[i] = tp->tree[tp->size -i -1]; }

	return tp;
}

void
del_linked_topol (linked_topol tp)
{
	int i;

	if (!tp) return;

	for (i=tp->size-1; i >= 0; i--) if (tp->tree[i]) del_topology (tp->tree[i]);
	if (tp->tree) free (tp->tree);

	free (tp);
}

void 
shuffle_topology (topology t, phylogeny *segment, linked_topol tp, split_space split)
{
//	int i, n_shuffle = biomc2_rng_uniform_int (random_number, 3);
	int i, n_shuffle = 0; 

	copy_topology_from_topology (tp->bkp[BKPrandom], t);
	for (i=0; i <= n_shuffle; i++) spr_apply_random_move_spr (tp->bkp[BKPrandom], true);

	copy_topology_from_topology_mapping (t, tp->bkp[BKPrandom], split);
	swap_likelihood (t, segment, tp->bkp[BKPrandom]);
}

void
pre_add_remove_topology (topology t, phylogeny *segment, linked_topol tp)
{
  clear_topology (t);
	copy_topology_from_topology (tp->bkp[BKPcopy], t);
  link_current_to_accepted (t, segment);
}

void
add_topology_forward (topology t, phylogeny *segment, linked_topol tp, split_space split)
{
  /* ( [t ] ) ==> ( [t->prev] [t ] ) segments 
	 * ( [T1] ) ==> ( [  T1*  ] [T1] ) topology */
	int i, mid = t->last_segment - t->first_segment;
  (void) split;

	insert_topology_left (t, tp);

	/* sort split point */
	if (mid > 1) {
		t->prev->first_segment = t->first_segment;
//		mid = (int) ((double)(mid)/2.);
//		i = mid + biomc2_rng_uniform_int (random_number, mid);
		i = biomc2_rng_uniform_int (random_number, mid);
		t->prev->last_segment = t->first_segment + i;
		t->first_segment += (i+1);
	}
	else {
		t->prev->first_segment = t->prev->last_segment = t->first_segment;
		t->first_segment = t->last_segment;
	}

	/* copy topology */ 
	copy_topology_from_topology (t->prev, t);
	clear_topology (t->prev);
	link_current_to_accepted (t->prev, segment);
//	shuffle_topology (t->prev, segment, tp, split);
//	spr_apply_random_move_spr (t->prev, false);
	(*branch_swap) (t->prev, false);
}

void
accept_add_topology_forward (topology t, phylogeny *segment)
{
	clear_topology (t->prev);
	link_accepted_to_current (t->prev, segment);
	sum_segment_likelihood (t, segment);
	t->likelihood_accepted = t->likelihood_current;
}

void
reject_add_topology_forward (topology t, phylogeny *segment, linked_topol tp)
{
  (void) segment;
//	unswap_likelihood (t->prev, segment, tp->bkp[BKPrandom]);
	t->first_segment = t->prev->first_segment;
	remove_topology (t->prev, tp);
}

void
add_topology_backward (topology t, phylogeny *segment, linked_topol tp, split_space split)
{
  /* ( [t ] ) ==> ( [t ] [t->next] ) 
	 * ( [T1] ) ==> ( [T1] [  T1*  ] ) */
	int i, mid = t->last_segment - t->first_segment;
  (void) split;

	insert_topology_right (t, tp);

	/* sort split point */
	if (mid > 1) {
		t->next->last_segment = t->last_segment;
//		mid = (int) ((double)(mid)/2.);
//		i = mid + biomc2_rng_uniform_int (random_number, mid);
		i = biomc2_rng_uniform_int (random_number, mid);
		t->next->first_segment = t->last_segment - i;
		t->last_segment -= (i+1);
	}
	else {
		t->next->first_segment = t->next->last_segment = t->last_segment;
		t->last_segment = t->first_segment;
	}

	/* copy topology */ 
	copy_topology_from_topology (t->next, t);
  clear_topology (t->next);
	link_current_to_accepted (t->next, segment);
//	shuffle_topology (t->next, segment, tp, split);
//	spr_apply_random_move_spr (t->next, false);
	(*branch_swap) (t->next, false);
}

void
accept_add_topology_backward (topology t, phylogeny *segment)
{
  clear_topology (t->next);
	link_accepted_to_current (t->next, segment);
	sum_segment_likelihood (t, segment);
	t->likelihood_accepted = t->likelihood_current;
}

void
reject_add_topology_backward (topology t, phylogeny *segment, linked_topol tp)
{
  (void) segment;
//	unswap_likelihood (t->next, segment, tp->bkp[BKPrandom]);
	t->last_segment = t->next->last_segment;
	remove_topology (t->next, tp);
}

void
remove_topology_forward (topology t, phylogeny *segment, linked_topol tp, split_space split) 
{
  (void) tp;
	/* ( [t->prev] ) ( [t] )  ==> ( [t->prev] ) ( [t->prev] ) */
	copy_topology_from_topology_mapping (t, t->prev, split);
	swap_likelihood (t, segment, t->prev);
}

void
accept_remove_topology_forward (topology t, phylogeny *segment, linked_topol tp)
{
	link_accepted_to_current (t, segment);
	clear_topology (t);
	t->likelihood_accepted += t->prev->likelihood_accepted;
	t->first_segment = t->prev->first_segment;
	remove_topology (t->prev, tp);
}

void
reject_remove_topology_forward (topology t, phylogeny *segment, linked_topol tp)
{
	unswap_likelihood (t, segment, t->prev);
	copy_topology_from_topology (t, tp->bkp[BKPcopy]);
	clear_topology (t);
}

void
remove_topology_backward (topology t, phylogeny *segment, linked_topol tp, split_space split) 
{
  (void) tp;
	/* ( [t] ) ( [t->next] ) ==> ( [t->next] ) ([t->next] ) */
	copy_topology_from_topology_mapping (t, t->next, split);
	swap_likelihood (t, segment, t->next);
}

void
accept_remove_topology_backward (topology t, phylogeny *segment, linked_topol tp)
{
	link_accepted_to_current (t, segment);
	clear_topology (t);
	t->likelihood_accepted += t->next->likelihood_accepted;
	t->last_segment = t->next->last_segment;
	remove_topology (t->next, tp);
}

void
reject_remove_topology_backward (topology t, phylogeny *segment, linked_topol tp)
{
	unswap_likelihood (t, segment, t->next);
	copy_topology_from_topology (t, tp->bkp[BKPcopy]);
	clear_topology (t);
}

void
shift_topology_forward (topology t, phylogeny *segment, linked_topol tp, split_space split)
{
  /* ( [t->prev] ) ( [t ] ) ==> ( [t->prev->prev] ) ( [t->prev] [t ] )  link structure
	 * ( [  T1   ] ) ( [T2] ) ==> ( [      T1     ] ) ( [  T1   ] [T2] )  topology     */
	int i;

	insert_topology_left (t, tp);

	/* sort split point */
	t->prev->first_segment = t->first_segment;
	i = biomc2_rng_uniform_int (random_number, (t->last_segment - t->first_segment));
	t->prev->last_segment = t->first_segment + i;
	t->first_segment += (i+1);
	
	/* copy topology preserving partial likelihood info (done/undone) */
	copy_topology_from_topology (t->prev, t);
	copy_topology_from_topology_mapping (t->prev, t->prev->prev, split);
	swap_likelihood (t->prev, segment, t->prev->prev);
	link_current_to_accepted (t->prev, segment);
}

void
accept_shift_topology_forward (topology t, phylogeny *segment, linked_topol tp)
{
	link_accepted_to_current (t->prev, segment);
	sum_segment_likelihood (t, segment);
	t->likelihood_accepted = t->likelihood_current;
	t->prev->prev->likelihood_accepted += t->prev->likelihood_accepted;
	t->prev->prev->last_segment = t->prev->last_segment;
	remove_topology (t->prev, tp);
}

void
reject_shift_topology_forward (topology t, phylogeny *segment, linked_topol tp)
{
	unswap_likelihood (t->prev, segment, t->prev->prev);
	t->first_segment = t->prev->first_segment;
	remove_topology (t->prev, tp);
}

void
shift_topology_backward (topology t, phylogeny *segment, linked_topol tp, split_space split)
{
  /* ( [t ] ) ( [t->next] ) ==> ( [t ] [t->next] ) ( [t->next->next] ) link structure
	 * ( [T1] ) ( [   T2  ] ) ==> ( [T1] [   T2  ] ) ( [      T2     ] ) topology     */
	int i;

	insert_topology_right (t, tp);

	/* sort split point */
	t->next->last_segment = t->last_segment;
	i = biomc2_rng_uniform_int (random_number, (t->last_segment - t->first_segment));
	t->next->first_segment = t->last_segment - i;
	t->last_segment -= (i+1);

	/* copy topology preserving partial likelihood info (done/undone) */
	copy_topology_from_topology (t->next, t);
	copy_topology_from_topology_mapping (t->next, t->next->next, split);
	swap_likelihood (t->next, segment, t->next->next);
	link_current_to_accepted (t->next, segment);
}

void
accept_shift_topology_backward (topology t, phylogeny *segment, linked_topol tp)
{
	link_accepted_to_current (t->next, segment);
	sum_segment_likelihood (t, segment);
	t->likelihood_accepted = t->likelihood_current;
	t->next->next->likelihood_accepted += t->next->likelihood_accepted;
	t->next->next->first_segment = t->next->first_segment;
	remove_topology (t->next, tp);
}

void
reject_shift_topology_backward (topology t, phylogeny *segment, linked_topol tp)
{
	unswap_likelihood (t->next, segment, t->next->next);
	t->last_segment = t->next->last_segment;
	remove_topology (t->next, tp);
}

void
swap_likelihood (topology old, phylogeny *segment, topology t)
{
	int i, j, k;

	/* t->nodelist[j]->id = j always;
	 * if ([j]->id == [j]->map_id) no swap is necessary;
	 * if ([j]->id == [ [j]->map_id ]->map_id) then swap will happen twice, 
	 * cancelling out, and so forth for larger circularities; that's why we need a
	 * mapping of likelihood vector beforehand (l_ptr[] is an extension of pivot to
	 * vectors) */

	for (j=t->nleaves; j < t->nnodes; j++) if (j != t->nodelist[j]->map_id) {
		k = j - t->nleaves;
		for (i=old->first_segment; i <= old->last_segment; i++) {
			segment[i]->l_ptr[k] = segment[i]->l[j];
		}
	}
	for (j=t->nleaves; j < t->nnodes; j++) if (j != t->nodelist[j]->map_id) {
		k = t->nodelist[j]->map_id - t->nleaves;
		for (i=old->first_segment; i <= old->last_segment; i++) {
			segment[i]->l[j] = segment[i]->l_ptr[k];
		}
	}
}

void
unswap_likelihood (topology old, phylogeny *segment, topology t)
{
	int i, j, k;

	/* revert swap_id and id */

	for (j=t->nleaves; j < t->nnodes; j++) if (j != t->nodelist[j]->map_id) {
		k = t->nodelist[j]->map_id - t->nleaves;
		for (i=old->first_segment; i <= old->last_segment; i++) {
			segment[i]->l_ptr[k] = segment[i]->l[j];
		}
	}
	for (j=t->nleaves; j < t->nnodes; j++) if (j != t->nodelist[j]->map_id) {
		k = j - t->nleaves;
		for (i=old->first_segment; i <= old->last_segment; i++) {
			segment[i]->l[j] = segment[i]->l_ptr[k];
		}
	}
}

void
link_current_to_accepted (topology t, phylogeny *segment)
{
	int i, j;
	t->likelihood_current = 0.; 
	for (i=t->first_segment; i <= t->last_segment; i++) { 
		t->likelihood_current += segment[i]->likelihood_accepted;
		segment[i]->likelihood_current = segment[i]->likelihood_accepted;
		for (j=t->nleaves; j < t->nnodes; j++) {
			segment[i]->l[j]->u_current = segment[i]->l[j]->u_accepted; 
			segment[i]->l[j]->d_current = segment[i]->l[j]->d_accepted; 
		}
	}
}

void
link_accepted_to_current (topology t, phylogeny *segment)
{
	int i, j;
	for (i=t->first_segment; i <= t->last_segment; i++) {
		segment[i]->likelihood_accepted = segment[i]->likelihood_current;
		for (j=t->nleaves; j < t->nnodes; j++) {
			segment[i]->l[j]->u_accepted = segment[i]->l[j]->u_current; 
			segment[i]->l[j]->d_accepted = segment[i]->l[j]->d_current; 
		}
	}
	t->likelihood_accepted = t->likelihood_current;
}

void
sum_segment_likelihood (topology t, phylogeny *segment)
{
	int i;
	t->likelihood_current = 0.; 
	for (i=t->first_segment; i <= t->last_segment; i++)
		t->likelihood_current += segment[i]->likelihood_accepted;
}

void
clear_topology (topology t)
{
  int i;
  for (i=t->nleaves; i < t->nnodes; i++) 
    t->nodelist[i]->u_done = t->nodelist[i]->d_done = true;
}
	
void
insert_topology_right (topology t, linked_topol tp)
{
	if (t->next) t->next->prev = tp->avail;
	else tp->last = tp->avail;
	tp->avail = tp->avail->next;
	tp->avail->prev->next = t->next;
	t->next = tp->avail->prev;
	t->next->prev = t;
	tp->avail->prev = NULL;
}

void
insert_topology_left (topology t, linked_topol tp)
{
	if (t->prev) t->prev->next = tp->avail;
	else tp->first = tp->avail;
	tp->avail = tp->avail->next; 
	tp->avail->prev->prev = t->prev; 
	t->prev = tp->avail->prev;
	t->prev->next = t;
	tp->avail->prev = NULL;
}

void
remove_topology (topology t, linked_topol tp)
{
	if (t->prev) t->prev->next = t->next;
	else tp->first = t->next;
	if (t->next) t->next->prev = t->prev;
	else tp->last = t->prev;

	tp->avail->prev = t;
	t->next = tp->avail;
	tp->avail = t;
	t->prev = NULL;
	/* DEBUG */
	t->root->left = NULL;
}

