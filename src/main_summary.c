#include "summary.h"

void print_mosaic_tree_frequencies (summary sm, mosaic m1, mosaic m2);
void print_recombinant_leaves (summary sm);
void print_centroid_mosaic_AISMpaper (summary sm, mosaic m);
void print_guidetree_index (summary sm, char *treefile, char *guidetreefile);

int
main (int argc, char **argv)
{
	clock_t time0, time1;
  summary sm;
  mosaic m1, m2;
  time0 = clock ();

  if ((argc != 3) && (argc != 4))
    biomc2_error ( " USAGE: %s <treefile> <distfile>\n OR\n %s <treefile> <distfile> <guide tree with true topols>", 
                   basename (argv[0]), basename (argv[0]));

  printf ("|\t %s\n| Recombination detection program based on segment-wise topology distance\n", ProgramVersion);
  printf ("|\t Leonardo de Oliveira Martins and Hirohisa Kishino\n");
  printf ("|\t\t This program is protected by the GPL license\n|\n");
  printf ("| Summary from posterior distribution\n\n");

	random_number = new_biomc2_rng ();

  sm = read_summary_from_file (argv[1], argv[2]);
	plot_post_recomb_freq (sm);
//  calculate_distance_matrix_bp_loc (sm);
  m1 = new_mosaic_from_centroid (sm);
  m2 = new_mosaic_from_quantile (sm, 0.5);

  print_mosaic_tree_frequencies (sm, m1, m2);
  print_recombinant_leaves (sm);
  //print_centroid_mosaic (sm, m);
  if (argc == 4) print_guidetree_index (sm, argv[1], argv[3]);
  
	printf ("original sample size: %d ; reduced: %d\n", sm->n_samples_original, sm->n_samples);

  del_mosaic (m1);
  del_mosaic (m2);
  del_summary (sm);
	time1 = clock ();
	fprintf (stdout, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	del_biomc2_rng (random_number);
	return EXIT_SUCCESS;
}

void
print_mosaic_tree_frequencies (summary sm, mosaic m1, mosaic m2)
{
  FILE *distfile = biomc2_fopen ("mosaicTreeFreq.txt", "w");
  int i, i1 = 0, i2 = 0, bp1, bp2, idx[sm->n_compact], tfreq[sm->n_compact], idx_second;
  double freq_second;
  char *s;
  empfreq e;

  for (i=0; i < sm->n_compact; i++) { idx[i]=i; tfreq[i]=0; }
  fprintf (distfile, "loc\t mapfreq\t maptre\t map2freq\t map2tre\t mapBP\t mosBP\t medBP\t recfreq\n");
  for (i=0; i < sm->n_segments-1; i++){
    bp1 = bp2 = 0;

    /* below we have bp++ since zero is reserved for _no_ recomb (bp is the tree ID) */
    if (i == m1->loc_idx[i1]) { bp1 = m1->tree[i1++]->i[0].idx; tfreq[bp1]++; bp1++; }
    if (i == m2->loc_idx[i2]) { bp2 = m2->tree[i2++]->i[0].idx; tfreq[bp2]++; bp2++; }
    tfreq[sm->tree_freq[i]->i[0].idx]++;

    /* it might be only one posterior topology */
    if (sm->tree_freq[i]->n > 1) {
      idx_second = sm->tree_freq[i]->i[1].idx + 1;
      freq_second = (double)(sm->tree_freq[i]->i[1].freq)/(double)(sm->n_samples_original);
      tfreq[sm->tree_freq[i]->i[1].idx]++;
    }
    else {
      idx_second = -1;
      freq_second = 0.;
    }

    fprintf (distfile, "%d\t %f\t %d\t %f\t %d\t %d\t %d\t %d\t %.6f\n", sm->breakpoint[i], 
             (double)(sm->tree_freq[i]->i[0].freq)/(double)(sm->n_samples_original),
             sm->tree_freq[i]->i[0].idx + 1, /* ID+1 whenever we _print_ tree ID */
             freq_second, idx_second, /* can be second most freq topol or ZERO */
             (sm->tree_freq[i]->i[0].idx == sm->tree_freq[i+1]->i[0].idx) ? (0) : (1),
             bp1, bp2,  /* index > 0 to avoid confusion with no-recomb */
             ((double) sm->recomb_freq[i]) / ((double) sm->n_samples_original));
  }

  tfreq[m1->tree[i1]->i[0].idx]++;
  tfreq[m2->tree[i2]->i[0].idx]++;
  tfreq[sm->tree_freq[i]->i[0].idx]++;
  if (sm->tree_freq[i]->n > 1) {
    idx_second = sm->tree_freq[i]->i[1].idx + 1;
    freq_second = (double)(sm->tree_freq[i]->i[1].freq)/(double)(sm->n_samples_original);
    tfreq[sm->tree_freq[i]->i[1].idx]++;
  }
  else {
    idx_second = -1;
    freq_second = 0.;
  }
  
  fprintf (distfile, "%d\t %f\t %d\t %f\t %d\t %d\t %d\t %d\t %.6f\n", sm->n_sites, 
           (double)(sm->tree_freq[i]->i[0].freq)/(double)(sm->n_samples_original),
           sm->tree_freq[i]->i[0].idx + 1, 
           freq_second, idx_second, 0, /* can be second most freq topol or ZERO ( last term is mapBP which is zero since there's no BP) */
           m1->tree[i1]->i[0].idx + 1, m2->tree[i2]->i[0].idx + 1, 0.);
  fclose (distfile);

  e = new_empfreq_from_int_weighted (idx, sm->n_compact, tfreq);
 
  distfile = biomc2_fopen ("mosaicTrees.tre", "w");
  fprintf (distfile, "#NEXUS\n\nBegin trees;\n Translate\n");
  fprintf (distfile, "\t1  %s", sm->taxlabel[0]);
  for (i=1; i < sm->t_compact[0]->nleaves; i++)
    fprintf (distfile, ",\n\t%d  %s", i+1, sm->taxlabel[i]);
  fprintf (distfile, "\n;\n");
  for (i=0; i < e->n; i++) {
    s = topology_to_string (sm->t_compact[ e->i[i].idx ]);
    fprintf (distfile, "tree t%d = %s\n",e->i[i].idx + 1, s);
    free (s);
  }
  fprintf (distfile, "\nEnd;\n");
  fclose (distfile);
  del_empfreq (e);
}

void
print_recombinant_leaves (summary sm)
{
  FILE *recfile = biomc2_fopen ("recSequences.txt", "w");
  int i, j, dist = 0, nleaves = sm->t_compact[0]->nleaves;
  int rec_sp[sm->t_compact[0]->nleaves], min = 1e6, max = 0;
  char *s;
  split_space split;
  (void) dist;

  split = new_split_space (nleaves);

  for (j=0; j < nleaves; j++) rec_sp[j] = 0;

  for (i=0; i < sm->n_segments-1; i++) if (sm->tree_freq[i]->i[0].idx != sm->tree_freq[i+1]->i[0].idx) {
    setup_split_store_recomb_info (split); /* extension to SPR algorithm that records pruned leaves */
    dist = dSPR_topology_l (sm->t_compact[sm->tree_freq[i  ]->i[0].idx], 
                            sm->t_compact[sm->tree_freq[i+1]->i[0].idx], split, &(split->inodes));

    fprintf (recfile, "BREAKPOINT "); /* for each breakpoint prints list of recomb leaves */
    for (j=0; j < nleaves; j++) if (split->leaf_recomb[j]) {
      fprintf (recfile, "%s\t", sm->taxlabel[j]); 
      rec_sp[j]++;
    }
    fprintf (recfile, "\n");
  }

  min = max = rec_sp[0];
  for (j=1; j < nleaves; j++) {
    if (rec_sp[j] < min) min = rec_sp[j];
    if (rec_sp[j] > max) max = rec_sp[j];
  }
  max -= min; if (!max) max = 1; /* avoid NaN */
  for (j=0; j < nleaves; j++) /* at the end, print frequency and standardized order of recombinant leaves */
    fprintf (recfile, "TOTAL %s\t%d\t%lf\n", sm->taxlabel[j], rec_sp[j], (double)(rec_sp[j] - min)/(double)(max));
  
  fclose (recfile);
  del_split_space (split);
 
  /* NEXUS tree file with MAP topologies only (for paup MAST calculation, e.g.) */
  recfile = biomc2_fopen ("mosaicTrees-MAP.tre", "w");
  fprintf (recfile, "#NEXUS\n\nBegin trees;\n Translate\n");
  fprintf (recfile, "\t1  %s", sm->taxlabel[0]);
  for (i=1; i < sm->t_compact[0]->nleaves; i++)
    fprintf (recfile, ",\n\t%d  %s", i+1, sm->taxlabel[i]);
  fprintf (recfile, "\n;\n");
  for (i=0; i < sm->n_segments-1; i++) if (sm->tree_freq[i]->i[0].idx != sm->tree_freq[i+1]->i[0].idx) {
    s = topology_to_string (sm->t_compact[ sm->tree_freq[i]->i[0].idx ]);
    fprintf (recfile, "tree t%d = %s\n", sm->tree_freq[i]->i[0].idx + 1, s);
    free (s);
  }
  s = topology_to_string (sm->t_compact[ sm->tree_freq[i]->i[0].idx ]); /* last segment */
  fprintf (recfile, "tree t%d = %s\n", sm->tree_freq[i]->i[0].idx + 1, s);
  free (s);
  fprintf (recfile, "\nEnd;\n");
  fclose (recfile);
}

void
print_centroid_mosaic_AISMpaper (summary sm, mosaic m)
{
  int i;
  printf ("OUTMOSAIC ");
  for (i=0; i < m->n; i++) {
    if (i) printf ("-%d", sm->breakpoint[m->loc_idx[i]]);
    else   printf ("%d", sm->breakpoint[m->loc_idx[i]]);
  }
  printf ("  ");
  for (i=0; i <= m->n; i++) {
    if (i) printf ("-%d", m->tree[i]->i[0].idx);
    else   printf ("%d", m->tree[i]->i[0].idx);
  }
  printf ("\n");
  printf ("OUTPOSTFREQ ");
  for (i=0; i < sm->n_segments-1; i++) {
    if (i) printf ("-%d", sm->breakpoint[i]);
    else   printf ("%d", sm->breakpoint[i]);
  }
  printf ("  ");
  for (i=0; i < sm->n_segments-1; i++) {
    if (i) printf ("-%6f", (double)(sm->recomb_freq[i])/(double)(sm->n_samples_original));
    else   printf ("%6f", (double)(sm->recomb_freq[i])/(double)(sm->n_samples_original));
  }
  printf ("\n");
}

void
print_guidetree_index (summary sm, char *treefile, char *guidetreefile)
{
  nexus_treespace nex, guidenex;
  int zero = 0, i, j, *label_order;
  bool found;
  topology *guidet;

  nex      = read_nexus_treespace_file (treefile);
  guidenex = read_nexus_treespace_file (guidetreefile);

  if (guidenex->nleaves != nex->nleaves) biomc2_error ( "number of taxa disagree between posterior and guide tree files");

  /* TODO leo @ 20090531: that's unnecessary to store them all in a vec, since I only lookup one at a time */
  guidet = (topology*) biomc2_malloc (guidenex->ntrees * sizeof (topology));

  /* order leaves (id numbers) based on nex->taxlabel_hash; mapping is stored in label_order */
  /*  OBS: here we index the leaves of guidetree based on the order in which they appear in treefile */
  /* TODO leo @ 20090605: These calls are innocuous, since for same treespace hashvector will have leaf IDs in order
   * (maybe some problem with different roots for distinct trees, wich order_nexus function fixes indeed)... */
  label_order = (int*) biomc2_malloc (guidenex->nleaves * sizeof (int));
  map_hashtable_order (label_order, nex->taxlabel_hash, guidenex->taxlabel, guidenex->nleaves);
  order_nexus_treespace_id (guidenex, label_order);

  printf ("OUTGUIDETREE ");
  for (i=0; i < guidenex->ntrees; i++) {
    guidet[i] = new_topology (guidenex->nleaves);
    copy_topology_from_nexus_tree (guidet[i], guidenex->T[i]);

    found = false;
    /* search for topology guidet in summary table of trees */
    for (j=0; (j < sm->n_compact) && (!found); j++) 
      if ((! dSPR_topology_l (sm->t_compact[j], guidet[i], sm->split, &zero))) {
        if (i) printf ("-%d",j);
        else   printf ("%d",j);
        found = true;
      }
    if (!found) { /* guide tree not found on posterior sample-> return number larger than any index */
      if (i) printf ("-%d",sm->n_compact);
      else   printf ("%d",sm->n_compact);
    }
  }

  printf ("\n");

  for (i=guidenex->ntrees-1; i >= 0; i--) if(guidet[i]) del_topology(guidet[i]);
  free (guidet);
  if (label_order) free (label_order);
  del_nexus_treespace (guidenex);
  del_nexus_treespace (nex);
}
 
