#include "run_sampler.h"

int
main (int argc, char **argv)
{
	clock_t time0, time1;
	
  time0 = clock ();

  if (argc != 2)
		biomc2_error ( " USAGE: %s <file>", basename (argv[0]));

  printf ("|\t %s\n| Recombination detection program based on segment-wise topology distance\n", ProgramVersion);
  printf ("|\t Leonardo de Oliveira Martins and Hirohisa Kishino\n|\n");
  printf ("| This program is protected by the GPL license\n\n");

	random_number = new_biomc2_rng ();

  run_sampler (argv[1]);
	
	time1 = clock ();
	fprintf (stderr, "timing: %.8f secs\n", (double)(time1-time0)/(double)CLOCKS_PER_SEC);

	del_biomc2_rng (random_number);
	return EXIT_SUCCESS;
}
