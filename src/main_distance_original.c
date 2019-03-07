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
	int j, count=0, spr1=0, no_replicates, recomb, overN = 0, underN = 0;
	bool success;
	double overest = 0., underest = 0., overvar = 0., undervar = 0., mean, var;
	time0 = clock ();
	
	if (argc != 3)
		biomc2_error (" USAGE: %s <tree file> <number of replicates>", basename (argv[0]));

	sscanf (argv[2], " %d ", &no_replicates);

	random_number = new_biomc2_rng ();
	
	treespace = read_nexus_treespace_file (argv[1]);
  j = biomc2_rng_uniform_int (random_number, treespace->ntrees);
  
  t1 = new_topology (treespace->nleaves);
	t2 = new_topology (treespace->nleaves);

  if (recomb > treespace->nleaves - 3) recomb = treespace->nleaves - 3;
	
	split = new_split_space (treespace->nleaves);

	copy_topology_from_nexus_tree (t1, treespace->T[j]);

//	s = topology_to_string_with_name_generate_length (t1, treespace->taxlabel, 0.2, 0.8);
//	s = topology_to_string (t1);
//	printf ("[ 0\t%d\tSPR 0 ] %s\n", recomb, s);
//	printf ("#NEXUS\nbegin trees;\n\ntree PAUP = [&U] %s\n", s);
//	free (s);

	time2 = 0;

  for (j = 0; j < no_replicates; j++) {
    recomb = biomc2_rng_uniform_int (random_number, treespace->ntrees);
    copy_topology_from_topology (t2, t1);
		success = spr_apply_multiple_moves (t1, recomb);

    if (success) {
			time1 = times (&times_buf);
			setup_split_store_recomb_info (split);
			spr1 = dSPR_topology_l (t2, t1, split, &(split->inodes));
			time2 += (times (&times_buf) - time1);

      
//			s = topology_to_string_with_name_generate_length (t1, treespace->taxlabel, 0.2, 0.8);
			s = topology_to_string (t1);
			printf ("[ %d\t%d\tSPR %d ] %s\n", j+1, recomb, spr1, s);
//      printf ("tree PAUP = [&U] %s\n", s);
			free (s);

//			printf ("%d\t%d\t%d\n", treespace->nleaves, recomb, spr1);

			if (recomb <= spr1) {
				overest += (spr1 - recomb);
				overvar += ((spr1 - recomb) * (spr1 - recomb));
				overN++;  
			}
			if (recomb >= spr1) {
				underest += (recomb - spr1);
				undervar += ((recomb - spr1) * (recomb - spr1));
				underN++;  
			}
			count++;
		}
	}

//	printf ("end;\n");

	fprintf (stderr, "%d\t%d\t", treespace->nleaves, recomb);

	mean = (overest - underest)/(double)(count);
	var  = ((overest - underest)*(overest - underest))/(double)(count);
	var  = overvar + undervar - var;
	if (count > 1) { var /= (double)(count - 1); var = sqrt (var); }
	else var = 0.;
	fprintf (stderr, "%.6f %.6f  ", mean, var);

	mean = (overest + underest)/(double)(count);
	var  = ((overest + underest)*(overest + underest))/(double)(count);
	var  = overvar + undervar - var;
	if (count > 1) { var /= (double)(count - 1); var = sqrt (var); }
	else var = 0.;
	fprintf (stderr, "%.6f %.6f  ", mean, var);

	if (overN) {
		mean = overest/(double)(overN);
		var  = (overest * overest)/(double)(overN);
		var  = overvar - var;
		if (overN > 1) { var /= (double)(overN - 1); var = sqrt (var); }
		else var = 0.;
		fprintf (stderr, "%.6f %.6f  ", mean, var);
	}
	else fprintf (stderr, "%.6f %.6f  ", 0., 0.);

	if (underN) {
		mean = underest/(double)(underN);
		var  = (underest * underest)/(double)(underN);
		var  = undervar - var;
		if (underN > 1) { var /= (double)(underN - 1); var = sqrt (var); }
		else var = 0.;
		fprintf (stderr, "%.6f %.6f  ", mean, var);
	}
	else fprintf (stderr, "%.6f %.6f  ", 0., 0.);
 

	del_topology (t1);
	del_topology (t2);
  
	del_nexus_treespace (treespace);
	del_split_space (split);
	del_biomc2_rng (random_number);

	time1 = clock ();

	fprintf (stderr, ". %.6fs ", (double)(time1-time0)/(double)CLOCKS_PER_SEC);
	fprintf (stderr, ". %.6fs ", (double)(time2)/(double) (sysconf (_SC_CLK_TCK)));
	fprintf (stderr, ". %d trees\n", count);
	
  return EXIT_SUCCESS;
}
