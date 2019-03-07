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

/*! \file hashtable.c
 *  \brief double hashing open-address hash table using strings as key
 *
 * Hash tables allow us to search for the position of a key (taxa name) without
 * scanning the whole vector (like in sequencial search).
 * This code is derived from the software DCM3, released under the GPL license.
 * The original copyright is 
 *
 * Copyright (C) 2004 The University of Texas at Austin
*/

#include "hashtable.h"

/*! \brief  hash function
 *
 * This is called before hash1() and hash2(). It sets the 
 * value of h to be used by hash1() and hash2() and saves excessive calls
 * to the hash function.
 */
unsigned long hash(hashtable ht, const char* str);

/*! \brief  hash function */
inline unsigned long hash1 (hashtable ht);

/*! \brief  hash function */
inline unsigned long hash2 (hashtable ht);

hashtable 
new_hashtable (int size)
{
	int i; 
	double x;
	hashtable ht;
  
	ht = (hashtable) biomc2_malloc (sizeof (struct hashtable_struct));

	/* setting hashtable to be size of power of 2
	 * this is required so that proper open addressing can occur,
	 * i.e. eventually every entry in the table is considered
	 */
	x = ceil (log (size) / log (2));
	size = pow (2, x+1);

  ht->size = size;
	
  /* allocate memory for table */
  ht->table = (hashtable_item*) biomc2_malloc(ht->size * sizeof (hashtable_item)); 

  /* initialize to NULLS*/
  for(i = 0; i < ht->size; i++) ht->table[i] = 0; 
  ht->P = 2147483647; /* initialize P (large prime)*/
  ht->probelength = 0;
  
  /*initialize hash1 and hash2 variables*/
  srand (time(0)); /* the GSL library would be overkill... */
  ht->a1 = rand() * (ht->P - 1) + 1;
  ht->a2 = rand() * (ht->P - 1) + 1;
  ht->b1 = rand() * ht->P;
  ht->b2 = rand() * ht->P;

  return ht;
}

void 
del_hashtable (hashtable ht)
{
	int i;
	for(i=ht->size-1; i>=0; i--)
		if (ht->table[i]) { 
			if(ht->table[i]->key) free (ht->table[i]->key); 
			free (ht->table[i]); 
		}
	if (ht->table) free (ht->table);
	if (ht) free (ht);
}


unsigned long 
hash (hashtable ht, const char* key) 
{
	unsigned long g;
	
	ht->h = 0;
	while (*key) {
		ht->h = (ht->h << 4) + *key++;
		g = ht->h & 0xF0000000L;
		if (g) ht->h ^= g >> 24;
		ht->h &= ~g;
  }
	return (ht->h % ht->P); /*P_ is a large prime*/
}

unsigned long 
hash1 (hashtable ht) 
{
	return ((((ht->a1 * ht->h) + ht->b1) % ht->P) % ht->size) ;
}

unsigned long 
hash2 (hashtable ht) 
{
  return ((((ht->a2 * ht->h + ht->b2) % ht->P) % (ht->size - 3)) | 1);
}

void 
insert_hashtable (hashtable ht, char* key, int value) 
{
	unsigned long i;
	int h1, h2;
	
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);

  ht->probelength = 0;
  
  for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
		ht->probelength++;
		if (!strcmp (ht->table[i]->key, key)) /* key was already inserted */
			return; 
	}
	/* alloc space for new key */
  ht->table[i] = biomc2_malloc (sizeof (struct hashtable_item_struct));
  ht->table[i]->key = (char*) biomc2_malloc((strlen (key)+1) * sizeof (char));
  strcpy(ht->table[i]->key, key);
  ht->table[i]->value = value;
  return;
}

int 
lookup_hashtable (hashtable ht, char* key) 
{
	unsigned long i;
	int h1, h2;
	
  hash (ht, key);
  h1 = hash1 (ht);
  h2 = hash2 (ht);
	
	ht->probelength = 0;

	for (i = h1; ht->table[i]; i = (i + h2) % ht->size) {
		ht->probelength++;
		if (!(ht->table[i])) return -1;
    else if ( (ht->table[i]) && (!strcmp (ht->table[i]->key, key)) ) 
      return ht->table[i]->value;
  }
  return -2;
}

int
binary_search_next_available (char **string_vector, int size)
{
	int first=0, last=size-1, middle;
	while (first < last) {
		middle = (int)((last+first)/2);
		if (string_vector[middle]) first = middle + 1;
		else last = middle;
	}
	if (!string_vector[last]) return last;
	else biomc2_error ("string_vector full; no space available\n");
  return 0;
}

