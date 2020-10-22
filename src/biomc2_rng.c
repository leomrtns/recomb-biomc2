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

/*! \file biomc2_rng.c
 * \brief Random number generation, from the Gnu Scientific Library (GSL)
 */

#include "biomc2_rng.h"

/*! \brief formula for Tausworthe number generation */
#define TAUSWORTHE(s,a,b,c,d) (((s &c) <<d) &RNG_MAX_VALUE) ^ ((((s <<a) &RNG_MAX_VALUE)^s) >>b)

/*! \brief returns a random number between 0 and RNG_MAX_VALUE */ 
inline unsigned long int biomc2_rng_get (biomc2_rng r);
/*! \brief Sets the new seed \f$s\f$ for the random generator. */
inline void biomc2_rng_set (biomc2_rng r, unsigned long int s);

double
biomc2_rng_uniform (biomc2_rng r)
{
  return (biomc2_rng_get (r) / 4294967296.0);
}

double
biomc2_rng_uniform_pos (biomc2_rng r)
{
	double x;
	do { x = (biomc2_rng_get (r) / 4294967296.0); } while (!x);
  return x;
}

unsigned long int
biomc2_rng_uniform_int (biomc2_rng r, unsigned long int n)
{
	unsigned long int scale = RNG_MAX_VALUE/n;
	unsigned long int k;

  do { k = biomc2_rng_get (r)/scale; } while (k >= n);
  return k;
}

inline unsigned long int
biomc2_rng_get (biomc2_rng r)
{
  r->s1 = TAUSWORTHE (r->s1, 13, 19, 4294967294UL, 12);
  r->s2 = TAUSWORTHE (r->s2, 2, 25, 4294967288UL, 4);
  r->s3 = TAUSWORTHE (r->s3, 3, 11, 4294967280UL, 17);
  return (r->s1 ^ r->s2 ^ r->s3);
}

inline void
biomc2_rng_set (biomc2_rng r, unsigned long int s)
{
	r->s1 = ((69069 * s) & RNG_MAX_VALUE);
	if (r->s1 < 2) r->s1 += 2UL;
	r->s2 = ((69069 * r->s1) & RNG_MAX_VALUE);
	if (r->s2 < 8) r->s2 += 8UL;
	r->s3 = ((69069 * r->s2) & RNG_MAX_VALUE);
	if (r->s3 < 16) r->s3 += 16UL;

  /* "warm it up" */
	biomc2_rng_get (r);
	biomc2_rng_get (r);
	biomc2_rng_get (r);
	biomc2_rng_get (r);
	biomc2_rng_get (r);
	biomc2_rng_get (r);
  return;
}

biomc2_rng 
new_biomc2_rng (void)
{
	biomc2_rng r = (biomc2_rng) biomc2_malloc (sizeof (struct biomc2_rng_struct));
	biomc2_rng_update (r);
  return r;
}

void
del_biomc2_rng (biomc2_rng r)
{
	if (r) free (r);
	return;
}

void
biomc2_rng_update (biomc2_rng r)
{
  struct timeval now;
  gettimeofday (&now, NULL);
  biomc2_rng_set (r, (unsigned long int)(now.tv_usec));
	return;
}

