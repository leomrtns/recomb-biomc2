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


#include "summary.h"

void
write_breakpoint_profile (summary sm)
{
  FILE *out = biomc2_fopen ("breakpoint-profile.tsv", "w");
  int mode = sm->nCOP_freq->i[0].idx;
  fprintf (out, "group\tposition\tweight\tbreak_weight\tdistance_sum\tmodal_count\talignment_length\tmetric\n");
  for (int group=0; group < 2; group++) {
    long long total = 0;
    for (int j=0; j < sm->n_samples; j++)
      if (!group || sm->nCOP[j] == mode) total += sm->w_sample[j];
    for (int i=0; i < sm->n_coldspot; i++)
      fprintf (out, "%s\t%d\t%lld\t0\t0\t%d\t%d\t%s\n", group ? "modal" : "all",
               sm->coldspot[i], total, mode, sm->n_sites, sm->metric);
    for (int i=0; i < sm->n_segments-1; i++) {
      long long breaks = 0, distance = 0;
      for (int j=0; j < sm->n_samples; j++) {
        if (group && sm->nCOP[j] != mode) continue;
        int w = sm->w_sample[j], d = sm->dSPR[i][j];
        if (d) breaks += w;
        distance += (long long)w * d;
      }
      fprintf (out, "%s\t%d\t%lld\t%lld\t%lld\t%d\t%d\t%s\n", group ? "modal" : "all",
               sm->breakpoint[i], total, breaks, distance, mode, sm->n_sites, sm->metric);
    }
  }
  fclose (out);
}
