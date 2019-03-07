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

/*! \file spr.h 
 *  \brief Header file for spr.c  
 */

#ifndef _biomc2_spr_h
#define _biomc2_spr_h

#include "topology.h"


/* Global function prototypes */

/*! \brief Apply a series of SPR operations avoiding previously visited topologies with same euler tour
 * composition. */
bool spr_apply_multiple_moves_eulertour (topology p, int n_moves);

/*! \brief Apply a series of SPR operations avoiding previously used rune or regraft
 * branches. */
bool spr_apply_multiple_moves (topology p, int n_moves);

/*! \brief Apply one branch swapping, which may be spr_apply_random_move_spr()
 * or spr_apply_random_move_nni() according to run_sampler.c */
void (*branch_swap) (topology p, bool done);

/*! \brief Apply SPR operation restricted to internal nodes (unused). */
void spr_apply_random_move_internal (topology p, bool done);

/*! \brief Apply one Subtree Prune-Regraft (SPR) operation. */
void spr_apply_random_move_spr (topology p, bool done);

/*! \brief Apply one Nearest-Neighbour Interchange (NNI) operation. */
void spr_apply_random_move_nni (topology p, bool done);

/*! \brief Apply one branch swapping operation at predetermined branches (when
 * reverting an SPR, for instance). */
void spr_apply_move_at_nodes (topology p, topol_node prune, topol_node regraft, bool not_phylo);

#endif
