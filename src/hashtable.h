/* Copyright (C) 2006	Leonardo de Oliveira Martins
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

/*! \file hashtable.h
 * \brief Header file for hashtable.c 
 */

#ifndef _biomc2_hashtable_h_ 
#define _biomc2_hashtable_h_ 

#include "biomc2_rng.h"


typedef struct hashtable_struct* hashtable;

typedef struct hashtable_item_struct* hashtable_item;

/*! \brief key/value pair for hash table */
struct hashtable_item_struct {
	/*! \brief String (vector of char). */
  char* key; 
	/*! \brief Integer (position in vector where hashtable_item_struct::key can be
	 * found) */
  int value;
};
/*! \brief Hash table (vector indexed by strings). */
struct hashtable_struct 
{ 
	/*! \brief Table size. */
	int size;
	/*! \brief Number of collisions before empty slot is found. */
	int probelength;
	/*! \brief Value set by hash(). Used in hash1() and hash2()
	 * to avoid calling hash() again. */
	unsigned long h;
	unsigned long a1, /*!< \brief Random values used in hash functions. */
								a2, /*!< \brief Random values used in hash functions. */
								b1, /*!< \brief Random values used in hash functions. */
								b2, /*!< \brief Random values used in hash functions. */
								P;  /*!< \brief Random values used in hash functions. */
	/*! \brief Vector with key/value pairs. */
	hashtable_item* table;
};

/*! \brief Insert key/value pair into hashtable. */
void insert_hashtable (hashtable ht, char* key, int value);
/*! \brief Return location (value) of corresponding key (string) or negative
 * value if not found. */
int  lookup_hashtable (hashtable ht, char* key);
/*! \brief Create new hashtable of size elements. */
hashtable new_hashtable (int size);
/*! \brief Free hashtable space. */
void del_hashtable (hashtable ht);
/*! \brief Binary search over string vector for next unused position.
 *
 *  When we do not have access to a counter, we need to do a sequential search
 *  over the vector of the form "while (not_there) { step++; }". A binary search
 *  allow us to decrease the number of steps. Notice that this is silly since an
 *  integer storing last value would do the trick (there is no reason why we
 *  do not have access to a counter...). 
 */
int binary_search_next_available (char **string_vector, int size);
#endif
