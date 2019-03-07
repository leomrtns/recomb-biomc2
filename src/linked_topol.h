/*! 
 * Copyright (C) 2006	Leonardo de Oliveira Martins
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

/*! \file linked_topol.h 
 *
 *  \brief doubly-linked list of topologies, with info about first and last
 *  segments. Topology updates are handled by these functions, without
 *  likelihood calculation.
 */

#ifndef _biomc2_linked_topol_h
#define _biomc2_linked_topol_h
#include "likelihood.h"

#define BKPSIZE 4
typedef enum
{
	BKPcopy = 0,       /*!< [1] original topology, before Al-Awadhi update    */
	BKPrandom,         /*!< multiple branch swaps on original topology in [1] */
	BKPshuffle_copy,   /*!< [2] original topology, independent from or outside Al-Awadhi update */
	BKPshuffle_random  /*!< multiple branch swaps on original topology in [2] */ 
} BackupTopology;

typedef struct linked_topol_struct* linked_topol;

struct linked_topol_struct
{
	/*! \brief doubly-linked list of topologies in use and available */
	topology *tree;
	/*! \brief list size */
	int size;
	/*! \brief pointer to first topology in use */
	topology first;
	/*! \brief pointer to last topology in use */
	topology last;
	/*! \brief pointer to first topology available */
	topology avail;
  /*! \brief pointer to topologies used for backup */
  topology bkp[BKPSIZE];
	/*! \brief number of recombinations and break-points */
	int nSPR, nCOP;
};

/* Global function prototypes */

linked_topol new_linked_topol (int nsegments, int nleaves);
void del_linked_topol (linked_topol tp);

void pre_add_remove_topology (topology t, phylogeny *segment, linked_topol tp);

void add_topology_forward  (topology t, phylogeny *segment, linked_topol tp, split_space split);
void add_topology_backward (topology t, phylogeny *segment, linked_topol tp, split_space split);
void accept_add_topology_forward  (topology t, phylogeny *segment);
void accept_add_topology_backward (topology t, phylogeny *segment);
void reject_add_topology_forward  (topology t, phylogeny *segment, linked_topol tp);
void reject_add_topology_backward (topology t, phylogeny *segment, linked_topol tp);

void remove_topology_forward  (topology t, phylogeny *segment, linked_topol tp, split_space split);
void remove_topology_backward (topology t, phylogeny *segment, linked_topol tp, split_space split);
void accept_remove_topology_forward  (topology t, phylogeny *segment, linked_topol tp);
void accept_remove_topology_backward (topology t, phylogeny *segment, linked_topol tp );
void reject_remove_topology_forward  (topology t, phylogeny *segment, linked_topol tp);
void reject_remove_topology_backward (topology t, phylogeny *segment, linked_topol tp);

void shift_topology_forward  (topology rt, phylogeny *segment, linked_topol tp, split_space split);
void shift_topology_backward (topology lt, phylogeny *segment, linked_topol tp, split_space space);
void accept_shift_topology_forward  (topology rt, phylogeny *segment, linked_topol tp);
void accept_shift_topology_backward (topology lt, phylogeny *segment, linked_topol tp);
void reject_shift_topology_forward  (topology rt, phylogeny *segment, linked_topol tp);
void reject_shift_topology_backward (topology lt, phylogeny *segment, linked_topol tp);

void   swap_likelihood (topology old, phylogeny *segment, topology t);
void unswap_likelihood (topology old, phylogeny *segment, topology t);
void link_current_to_accepted (topology t, phylogeny *segment);
void link_accepted_to_current (topology t, phylogeny *segment);

void clear_topology (topology t);
void sum_segment_likelihood (topology t, phylogeny *segment);
void remove_topology (topology t, linked_topol tp);

#endif
