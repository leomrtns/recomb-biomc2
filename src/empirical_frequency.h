/*!
 * Copyright (C) 2006 Leonardo de Oliveira Martins
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

/*! \file empirical_frequency.h
 *  \brief Header file for empirical_frequency.c
 */

#ifndef _empirical_frequency_h
#define _empirical_frequency_h

#include "topology.h"

typedef struct empfreq_struct* empfreq;

typedef struct
{
	int freq;
	int idx;
} empfreq_element;

struct empfreq_struct
{
	empfreq_element *i;
	int n;
  /*! \brief Min value for index. */
  int min; 
  /*! \brief Max value for index. */
  int max;
};

void sort_empfreq_decreasing (empfreq ef);
void sort_empfreq_increasing (empfreq ef);

empfreq new_empfreq (int n_elements);
void del_empfreq (empfreq ef);

empfreq new_empfreq_from_int_weighted (int *vec, int n, int *weight);
empfreq new_empfreq_from_int (int *vec, int n);
int find_mode_int_weighted (int *vec, int n, int *weight);
int find_mode_int (int *vec, int n);


#endif
