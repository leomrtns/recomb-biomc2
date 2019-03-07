#include "biomc2.h"
#include "topology.h"


int remove_duplicate_trees (topology *p, int ntrees, split_space split);
void summarise_treespace (topology *p1, topology *p2, int distinct1, int distinct2, split_space split, nexus_treespace t1);

int
main (int argc, char **argv)
{
  clock_t time0, time1;
  nexus_treespace t1, t2;
  topology *p1, *p2;
  int j, ntrees1, ntrees2, distinct1, distinct2, *label_order;
  split_space split;

  time0 = clock ();

  if (argc != 3){
    biomc2_error ( " USAGE: %s <tree file 1> <tree file 2>", basename (argv[0]));
  }

  fprintf (stderr, "|\t %s\n| Recombination detection program based on segment-wise topology distance\n", ProgramVersion);
  fprintf (stderr, "|\t Leonardo de Oliveira Martins and Hirohisa Kishino\n|\n");
  fprintf (stderr, "| This program is protected by the GPL license \n|\n");
  fprintf (stderr, "|\t \t Distribution of SPR distances \n\n");

  /* read topologies information */
  t1 = read_nexus_treespace_file (argv[1]);
  p1 = (topology*) biomc2_malloc (t1->ntrees * sizeof (topology));
  split = new_split_space (t1->nleaves);

  label_order = (int*) biomc2_malloc (t1->nleaves * sizeof (int));
  /* order leaves (id numbers) based on t1->taxlabel_hash; mapping is stored in label_order */
  map_hashtable_order (label_order, t1->taxlabel_hash, t1->taxlabel, t1->nleaves);
  order_nexus_treespace_id (t1, label_order);

  for (j=0; j < t1->ntrees; j++) {
    p1[j] = new_topology (t1->nleaves);
    copy_topology_from_nexus_tree (p1[j], t1->T[j]);
  }
  ntrees1 = t1->ntrees;
  distinct1 = remove_duplicate_trees (p1, ntrees1, split);

  t2 = read_nexus_treespace_file (argv[2]);
  if (t2->nleaves != t1->nleaves) biomc2_error ( "number of taxa disagree between files ");
  p2 = (topology*) biomc2_malloc (t2->ntrees * sizeof (topology));
  /* order leaves (id numbers) based on t1->taxlabel_hash; mapping is stored in label_order */
  map_hashtable_order (label_order, t1->taxlabel_hash, t2->taxlabel, t1->nleaves);
  order_nexus_treespace_id (t2, label_order);

  for (j=0; j < t2->ntrees; j++) {
    p2[j] = new_topology (t2->nleaves);
    copy_topology_from_nexus_tree (p2[j], t2->T[j]);
  }
  ntrees2 = t2->ntrees;
  distinct2 = remove_duplicate_trees (p2, ntrees2, split);

  summarise_treespace (p1, p2, distinct1, distinct2, split, t1);

  for (j=ntrees1-1; j >= 0; j--) if(p1[j]) del_topology(p1[j]);
  free (p1);
  for (j=ntrees2-1; j >= 0; j--) if(p2[j]) del_topology(p2[j]);
  free (p2);
  if (label_order) free (label_order);
  del_split_space (split);
  del_nexus_treespace (t1);
  del_nexus_treespace (t2);

  time1 = clock ();
  fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

  return EXIT_SUCCESS;
}

int		
remove_duplicate_trees (topology *p, int ntrees, split_space split)
{
  int distinct = ntrees;
  int i, j, zero = 0;
  topology pivot;


  for (i=0; i < ntrees; i++)  p[i]->tfreq = 1;

  for (i=0; i < distinct-1; i++) {
    for (j=i+1; j < distinct; j++) 
      if ((! dSPR_topology_l (p[i], p[j], split, &zero))  ) {
        distinct--;
        pivot      = p[distinct];
        p[distinct] = p[j]; 
        p[j]       = pivot;
        p[i]->tfreq++;
        j--;
      }
  }

  return distinct;
}

void
summarise_treespace (topology *p1, topology *p2, int distinct1, int distinct2, split_space split, nexus_treespace t1)
{
	int histogram[p1[0]->nleaves], i, j, k, dist;
	int recomb_sp[p1[0]->nleaves], n_recomb = 0, freq;
	FILE *outfile;
	for (i=0; i < p1[0]->nleaves; i++) { histogram[i] = 0; recomb_sp[i] = 0; }

  for (i=0; i < distinct1; i++)
    for (j=0; j < distinct2; j++) {
			setup_split_store_recomb_info (split);
      dist = dSPR_topology_l (p1[i], p2[j], split, &(split->inodes));
			freq = (p1[i]->tfreq * p2[j]->tfreq);
			histogram[dist] += freq;
			if (dist) {
				n_recomb += freq;
				for (k=0; k < p1[0]->nleaves; k++) if (split->leaf_recomb[k])
					recomb_sp[k] += freq;
			}
		}

  for (i=0; i < p1[0]->nleaves; i++) 
    printf ("%d\t %d\n", i, histogram[i]);

	outfile = biomc2_fopen ("recomb_leaves.txt", "w");
	fprintf (outfile, "taxon \t rfreq\n");
	for (i=0; i < p1[0]->nleaves; i++) 
		fprintf (outfile, "%s \t %f\n", t1->taxlabel[i], ((double) recomb_sp[i])/((double) n_recomb) );
	fclose (outfile);
}

