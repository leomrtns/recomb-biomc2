# from http://www.gnu.org/software/hello/manual/automake/
# foo_CFLAGS is per-target; AM_CFLAGS is for all targets except the per-targets; 
# CFLAGS is user-defined outside configure, applied after AM_CFLAGS (or foo_CFLAGS)
biomc2_summarise_CPPFLAGS = $(DEPS_CFLAGS)
# LDADD is not user-defined, is a makefile-wide
biomc2_summarise_LDADD = $(DEPS_LIBS)
#AM_CFLAGS = $(cairo_CFLAGS)
#LIBS = $(cairo_LIBS)

common_headers = biomc2.h biomc2_rng.h hashtable.h nexus.h topology.h
common_src = 	biomc2.c biomc2_rng.c hashtable.c nexus.c nexus_alignment.c nexus_tree.c \
							topology.c topology_distance.c topology_split_partition.c

bmc2_headers = 	spr.h phylogeny.h likelihood.h linked_topol.h chain_data.h \
								update_chain.h convergence.h run_sampler.h
bmc2_src =	spr.c phylogeny.c likelihood.c linked_topol.c chain_data.c \
						update_chain_params.c update_chain_topology.c convergence.c \
						run_sampler.c

sumry_headers = summary.h empirical_frequency.h
sumry_src = summary.c empirical_frequency.c main_summary.c summary_plot.c summary_cluster.c

bin_PROGRAMS = biomc2 biomc2.summarise biomc2.distance biomc2.sprdist
biomc2_SOURCES = $(common_src) $(bmc2_src) $(common_headers) $(bmc2_headers) main.c
biomc2_summarise_SOURCES = $(common_src) $(sumry_src) $(common_headers) $(sumry_headers)
biomc2_sprdist_SOURCES = $(common_src) $(common_headers) main_sprdist.c
biomc2_distance_SOURCES = $(common_src) spr.c $(common_headers) spr.h main_distance.c

