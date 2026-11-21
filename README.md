# biomc2
The software **biomc2** is a hierarchical Bayesian procedure to detect recombination in DNA sequence alignments based on the incongruence between topologies caused by recombinational events.
Trees can span multiple 64-bit words. All segments use the same set of labelled tips.

The model is described in the following publications:
  * Leonardo de Oliveira Martins, Élcio Leal and Hirohisa Kishino (2008) Phylogenetic Detection of Recombination with a Bayesian Prior on the Distance between Trees. _PLoS ONE_ 3(7): e2651. [doi:10.1371/journal.pone.0002651](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0002651)

  * Leonardo de Oliveira Martins and Hirohisa Kishino (2010) Distribution of distances between topologies and its effect on detection of phylogenetic recombination. _Annals of the Institute of Statistical Mathematics_ 62(1): 145--159. [doi:10.1007/s10463-009-0259-8](http://www.springerlink.com/content/y550n510j810w8k2/)

## Installation
From a Git checkout, generate the build files with Autoconf and Automake. A source tarball produced by `make dist` already includes `configure`.

```sh
git clone git@github.com:leomrtns/recomb-biomc2.git
cd recomb-biomc2
autoreconf -fi
mkdir build && cd build
../configure --prefix="$HOME/.local"
make -j4
make check
make install
```

The four executables require a C99-capable compiler and the system maths library. Set `CC`, `CFLAGS` and `LDFLAGS` before configuring to select a toolchain. Python 3 is used by `make check`. Base R is optional and is used by `biomc2.plot`. For a staged installation use `make install DESTDIR="$staging_directory"`. A local Conda recipe is provided in `conda/`. From the repository root run `conda build -c conda-forge conda` in an environment with conda-build.

## Programs

This suite is composed of four programs: "biomc2", "biomc2.summarise", "biomc.sprdist" and "biomc2.distance", which are explained in more detail below. Despite my efforts in complying with the GNU standards, the usage of the programs is quite "unusual": the parameters are given in a specific order. If you simply run the program without any parameters, it will give a minimalistic  hint about its usage.
### biomc2
Usage: biomc2 control_file
This is the main program, that given a DNA alignment and a tree file will try to estimate the recombination break-points and the number of recombinations for each break-point. The alignment file should be in INTERLEAVED nexus format. Or at least have the word "INTERLEAVE", since as far as I understand the non-interleaved format is a special case (the whole sequence in one line).
        
The tree file is just an initial guess and should be in nexus format. If you feel unconfortable about the arbitrariness of this initial state you can give a tree file with lots of trees. The program will choose one randomly. All trees should be bifurcating (no politomies), with exception of the root node - since  unrooted trees are represented by a trifurcation, in parenthetic format. The  program can read rooted and unrooted trees, since in the end it will remove an eventual root node.

The only input parameter to the program is the control file name. In this file there will be information about the DNA and tree files, along with other model parameters. A template file can be found at doc/ctrl.biomc2 and an example is available at example/ctrl.biomc2 . 

Set `metric = dspr` (default), `rf` or `hdist` in the control file. The sampler applies the selected distance between neighbouring segment trees in its prior. Approximate dSPR uses the split-pruning algorithm adapted from genefam-dist. RF counts splits present in only one of the two trees (symmetric difference). Hdist sums the minimum assignment cost between splits, with the cost for a pair equal to the smaller Hamming distance across the two orientations. RF avoids assignment. Hdist solves one assignment on the unmatched edges. Approximate dSPR may require several rounds of pruning.

`maxdist` is the largest allowed neighbour distance in units of the selected metric. For example, an NNI gives RF = 2, so `maxdist = 1` with RF excludes all topology changes. An explicit zero selects the historical binary same/different prior. Omitting the cutoff uses the number of tips minus three. The choice and scale of distance affect inference. In particular, RF and Hdist have no interpretation as counts of recombination events.

The sampler writes `post.tree` and `post.dist` (or `prior.tree` and `prior.dist` for a prior run). The trace header records the metric. Each sample gives its iteration, log score, total tree distance, number of breakpoints, prior average rate and prior kappa, followed by the breakpoint indices and their distances. Each breakpoint index identifies a boundary in the second header line.
### biomc2.summarise
Usage: biomc2.summarise post.tree post.dist
This program reads the treefile and distfile from the MCMC sampler and summarises the posterior distribution of break-points. It changed a lot from the previous version, and is focused on determining a likely mosaic structure. The caveat is that it can no  longer summarise the posterior distribution of average rates per segment - but this was  not our main purpose in the first place. If necessary these and other values can be output with a little programming (on function sample_to_output_file() of file run_sampler.c).

The program uses two main approaches to estimate a break-point mosaic structure: the piece-wise median and the centroid sample. The piece-wise median is based on the idea that the posterior distribution of break-points will have several modes, one for each break-point, along the alignment. Thus we partition the distribution based on the CDF where the number of partitions equals the modal number of break-points, and for each partition we find its median value. The centroid sample is in fact the sample whose mosaic structure minimizes the total distance to all other samples. This distance algorithm will be published in the "Annals of the Institute of Statistical Mathematics". For each estimated mosaic structure we infer the topology as the MAP topology over all segments composing the non-recombinant regions.

The program writes the mosaic tables and trees, candidate recombinant leaves in `recSequences.txt`, and a plotting table in `breakpoint-profile.tsv`. Pruned leaf labels are inferred by dSPR between neighbouring MAP trees for all three sampler metrics. The file "mosaicTrees.tre" contains a list with all relevant topologies, named as character "t" followed by a numerical ID (used by mosaicTreeFreq.txt). The file "mosaicTreeFreq.txt" contains a table with 9 columns, where each row represents one segment:

  - loc: location of last site of segment;
  - mapfreq: posterior frequency of MAP (Maximum A Posteriori, the mode) topology for segment;
  - maptre: ID of MAP tree (according to mosaicTrees.tre);
  - map2freq: posterior frequency of second most frequent topology;
  - map2tre: tree ID of second most frequent topology;
  - mapBP: indicator of whether the MAP topologies changed or not between the segment and next;
  - mosBP: mosaic structure according to centroid sample (zero means no
    break-point and number other than zero is the tree ID up to this segment);
  - medBP: mosaic structure based on median (the values different from zero indicate the estimated topology ID for segments up to this one);

  - recfreq: posterior probability of a breakpoint at this boundary.

Run `biomc2.plot` after `biomc2.summarise` to write `recomb_freq.pdf`. Optional positional arguments select the input TSV and output PDF. Each panel shows the posterior mean tree distance and breakpoint probability by alignment position. Orange bars give the central 95% interval within each equal-mass partition of the breakpoint distribution, with a diamond at its median. Rug marks indicate boundaries with zero observed breakpoints in that panel.

The table contains one row per candidate boundary for each group (`all` and draws with the modal breakpoint count, `modal`). Divide `distance_sum` or `break_weight` by `weight` to obtain the posterior mean distance or breakpoint probability. All three columns count the multiplicity of compacted samples. Historical traces without a metric field are read as dSPR.

If you are curious: if a third argument is given to biomc2.summarise with a tree file it will compare the trees in this file with the posterior trees and output the corresponding index (tree ID) in the posterior samples (I used this obscure option when comparing with simulations, where we know the true tree and want to see if the sampler can find it). Don't worry bout this option.


### biomc2.distance 
`biomc2.distance tree_file number_of_replicates` applies random SPR moves and reports the number applied, approximate dSPR and half the RF distance. Simulated trees are written to `spr.trees`.
### biomc2.sprdist
`biomc2.sprdist tree1 tree2 [dspr|rf|hdist]` compares the two sets of trees and prints a histogram weighted by their observed frequencies. With the default dSPR it also writes candidate recombinant leaf frequencies to `recomb_leaves.txt`.
## License
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2008-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)
