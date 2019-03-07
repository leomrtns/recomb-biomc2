#include "spr.h"

int
main (int argc, char **argv)
{
	clock_t time0, time1, time2;
	struct tms times_buf;
	nexus_treespace treespace;
	split_space split;
	topology t1, t2;
	char *s;
  FILE *stream;
	int j, count=0, spr1=0, no_replicates, recomb;
	bool success;
	time0 = clock ();
	
	if (argc != 3)
		biomc2_error (" USAGE: %s <tree file> <number of replicates>", basename (argv[0]));

	sscanf (argv[2], " %d ", &no_replicates);

	random_number = new_biomc2_rng ();
	
	treespace = read_nexus_treespace_file (argv[1]);
  j = biomc2_rng_uniform_int (random_number, treespace->ntrees);
  
  t1 = new_topology (treespace->nleaves);
	t2 = new_topology (treespace->nleaves);
	split = new_split_space (treespace->nleaves);

	copy_topology_from_nexus_tree (t1, treespace->T[j]);

  stream = fopen ("spr.trees", "w");
	s = topology_to_string (t1);
  fprintf (stream, "#NEXUS\nbegin trees;\n\ntree PAUP0 = [&U] %s\n", s);
	free (s);

	time2 = 0;

  for (j = 0; j < no_replicates; j++) {
    recomb = biomc2_rng_uniform_int (random_number, treespace->nleaves - 3) + 1;
    copy_topology_from_topology (t2, t1);
		success = spr_apply_multiple_moves (t1, recomb);

    if (success) {
			time1 = times (&times_buf);
			setup_split_store_recomb_info (split);
			spr1 = dSPR_topology_l (t2, t1, split, &(split->inodes));
			time2 += (times (&times_buf) - time1);
			printf ("%8d %8d %8d\n", recomb, spr1, split->rfdistance);

			s = topology_to_string (t1);
      fprintf (stream, "tree PAUP%d = [&U] %s\n", j+1, s);
			free (s);
		}
	}

	del_topology (t1);
	del_topology (t2);
  
  fprintf (stream, "end;\n");
  fclose (stream);

	del_nexus_treespace (treespace);
	del_split_space (split);
	del_biomc2_rng (random_number);

	time1 = clock ();

	fprintf (stderr, ". %.6fs ", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
	fprintf (stderr, ". %.6fs ", (double)(time2)/(double) (sysconf (_SC_CLK_TCK)));
	fprintf (stderr, ". %d trees\n", count);
	
  return EXIT_SUCCESS;
}
