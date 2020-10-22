# biomc2
The software **biomc2** is a hierarchical Bayesian procedure to detect recombination in DNA sequence alignments based on the incongruence between topologies caused by recombinational events.

The model is described in the following publications:
  * Leonardo de Oliveira Martins, Élcio Leal and Hirohisa Kishino (2008) Phylogenetic Detection of Recombination with a Bayesian Prior on the Distance between Trees. _PLoS ONE_ 3(7): e2651. [doi:10.1371/journal.pone.0002651](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0002651)

  * Leonardo de Oliveira Martins and Hirohisa Kishino (2010) Distribution of distances between topologies and its effect on detection of phylogenetic recombination. _Annals of the Institute of Statistical Mathematics_ 62(1): 145--159. [doi:10.1007/s10463-009-0259-8](http://www.springerlink.com/content/y550n510j810w8k2/)

For more information please check out http://www.biomcmc.org

## Installation

Some general instructions on how to install the package are described in the INSTALL file. 
The configuration is handled by autotools, so if you have trouble consider running `autoreconf` in the root directory (see below).
I've tested only on GNU/Debian linux machines (intel and AMD64), and MacOSX using fink (http://fink.sourceforge.net/) but it may work on Windows using cygwin 
(http://www.cygwin.com/). 

### quick installation guidelines

```[bash]
/home/simpson/$ git clone https://github.com/leomrtns/recomb-biomc2.git
/home/simpson/$ cd recomb-biomc2 && autoreconf
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../recomb-biomc2/configure  --enable-static-binary --prefix=${HOME}/local # or a place of your choice for installation
/home/simpson/$ make; make install
```

* The "build" directory is where the compilation will be done, and it's good practice to build the program in a separate directory. 
* The option "--enable-static-binary" will link the compiled files with "-static", so that you can copy the binaries to other machines with the same architecture/OS but which lack the libraries. 
* Another option specific to these programs is "--enable-optimization", which will include optimization flags when compiling. The program may become faster, but I don't know if it works for other architectures or older compilers. 
* An important parameter is "--prefix=SOMEWHERE", since the compiled programs will go to "SOMEWHERE/bin/" and so forth. 
The default is '/usr/local', but it may be better to install it locally. In the example above, it will go to `${HOME}/local`, where `${HOME}` is your home directory (`/home/simpson` in the example).
* In some cases, to compile it statically, you may need to configure with a command like (at least the `pkg_config` and `-static` are obligatory...):
```bash
PKG_CONFIG="pkg-config --static" CFLAGS="-pthread -static -flto -ffat-lto-objects -fPIC" ../biomc2-1.11/configure --enable-static-binary --enable-optimization
```
 * If `configure` complains about a missing library, you'll need to install them before running `configure` again.
In particular it depends on the GTK development library, for plotting the figure when running biomc2.summarise. 
 * If there is no `configure` file at all in the distribution, or it doesn't run, then you you need to install the `autotools` and rerun the configuration.
Both cases are shown below, if you can install them system-wide:

```[bash]
## 'bootstrap' the configuration files (needed when cloning from github):
/home/simpson/$ apt-get install autotools-dev autoconf automake libtool
/home/simpson/$ (cd recomb-biomc2 && autoreconf)  ## the parentheses avoid entering the directory afterwards
### install missing libraries
/home/simpson/$ apt-get install apt-get install libgtk2.0-dev libcairo2-dev
 ```
  * The libraries rely on `pkg-config` to find their location: if your `pkg-config` was installed through conda then you'd
better install the above libs via conda as well (or, you know, updating environmental variables etc)

### recommended software/libraries
* The R software for statistical computing (http://www.r-project.org/). It is useful for analysing the output of the program, plotting, etc. The following libraries for R are very handy:
  * ape (http://cran.r-project.org/src/contrib/Descriptions/ape.html), for plotting the trees, calculating consensus trees, even estimating population parameters;
  * coda (http://cran.r-project.org/src/contrib/Descriptions/coda.html), for analysing and diagnosing MCMC simulations;
  * bioconductor (http://www.bioconductor.org/), for analysis of genomic data in general (I use it for creating nice figures).
* MrBayes (http://mrbayes.csit.fsu.edu/), for analysing the posterior distribution of topologies. In the first version of our program, MrBayes was necessary. From version 1.9 on the program biomc2.summarise can calculate the MAP topologies.
* paml (http://abacus.gene.ucl.ac.uk/software/paml.html), for simulating alignments along a given phylogeny.
- paup (http://paup.csit.fsu.edu/), for generating random topologies. This is not free software, so you may try the 'rtree' function from the "ape" package in R. The program 'biomc2.distance' also generates trees, but they are correlated by the recombination distance.

## Programs

This suite is composed of four programs: "biomc2", "biomc2.summarise", "biomc.sprdist" and "biomc2.distance", which are explained in more detail below. Despite my efforts in complying with the GNU standards, the usage of the programs is quite "unusual": the parameters are given in a specific order. If you simply run the program without any parameters, it will give a minimalistic  hint about its usage.
### biomc2
Usage: biomc2 control_file
This is the main program, that given a DNA alignment and a tree file will try to estimate the recombination break-points and the number of recombinations for each break-point. The alignment file should be in INTERLEAVED nexus format. Or at least have the word "INTERLEAVE", since as far as I understand the non-interleaved format is a special case (the whole sequence in one line).
        
The tree file is just an initial guess and should be in nexus format. If you feel unconfortable about the arbitrariness of this initial state you can give a tree file with lots of trees. The program will choose one randomly. All trees should be bifurcating (no politomies), with exception of the root node - since  unrooted trees are represented by a trifurcation, in parenthetic format. The  program can read rooted and unrooted trees, since in the end it will remove an eventual root node.

The only input parameter to the program is the control file name. In this file there will be information about the DNA and tree files, along with other model parameters. A template file can be found at doc/ctrl.biomc2 and an example is available at example/ctrl.biomc2 . 

The program outputs to the screen progress information, consisting of the acceptance rates for each move and a point estimate of the number of recombinations nSPR and the number of break-points nCOP in the form [nSPR nCOP] . It also outputs the sampled values of each parameter to the files "post.tre" and "post.dist". If you are recovering the prior distribution (in other words, if you  set "prior=1" in control_file) then the output files will be named "prior.tre" and "prior.dist". (in the original version the program would produce way too many output files, restricting the maximum number of segments). 

The file "post.dist" has a variable number of columns and should be processed by program "biomc2.summarise" but the four first columns (after the third row) are respectively "iteration", "posterior probability" (proportional to; in log scale), "number of SPRs" (=lower bound on number of recomb events), "number of break-points" (one break-point can harbor more than one recomb event), "prior average rate" and "prior kappa of HKY". The file "post.tree" is a standard nexus tree file with trees mapped to segments and samples (you need "biomc2.summarise" to make sense out of them, though).


### biomc2.summarise
Usage: biomc2.summarise post.tre post.dist 
This program reads the treefile and distfile from the MCMC sampler and summarises the posterior distribution of break-points. It changed a lot from the previous version, and is focused on determining a likely mosaic structure. The caveat is that it can no  longer summarise the posterior distribution of average rates per segment - but this was  not our main purpose in the first place. If necessary these and other values can be output with a little programming (on function sample_to_output_file() of file run_sampler.c).

The program uses two main approaches to estimate a break-point mosaic structure: the piece-wise median and the centroid sample. The piece-wise median is based on the idea that the posterior distribution of break-points will have several modes, one for each break-point, along the alignment. Thus we partition the distribution based on the CDF where the number of partitions equals the median number of break-points, and for each partition we find its median value. The centroid sample is in fact the sample whose mosaic structure minimizes the total distance to all other samples. This distance algorithm will be published in the "Annals of the Institute of Statistical Mathematics". For each estimated mosaic structure we infer the topology as the MAP topology over all segments composing the non-recombinant regions. 

The program output is composed of three files: "mosaicTreeFreq.txt", "mosaicTrees.tre" and "recomb_freq.pdf". The file "mosaicTrees.tre" contains a list with all relevant topologies, named as character "t" followed by a numerical ID (used by mosaicTreeFreq.txt). The file "mosaicTreeFreq.txt" contains a table with 8 columns, where each row represents one segment:

  - loc: location of last site of segment;
  - mapfreq: posterior frequency of MAP (Maximum A Posteriori, the mode) topology for segment;
  - maptre: ID of MAP tree (according to mosaicTrees.tre);
  - map2freq: posterior frequency of second most frequent topology;
  - map2tre: tree ID of second most frequent topology;
  - mapBP: indicator of whether the MAP topologies changed or not between the segment and next;
  - mosBP: mosaic structure according to centroid sample (zero means no
    break-point and number other than zero is the tree ID up to this segment);
  - medBP: mosaic structure based on median (the values different from zero indicate the estimated topology ID for segments up to this one);

The file "recomb_freq.pdf" is a figure with the posterior distribution of (approx.) recombinations, where the upper panel uses all samples and lower panel  uses only the fraction of samples that exhibit the modal (most frequent) number of break-points. For each panel, we have at the top the mosaic based on the piece-wise median (together with 95% inter-quantile ranges). If this estimate does  not coincide with the "mosaicTreeFreq.txt" please use only the mosaicTreeFreq.txt one and report to me. I haven't verified the plotting functions... At the bottom we have blue and red bars: the blue bars are segments where in no sample we observed a break-point (cold spots). The red bars represent the 95% credible set (together, they are responsible for 95% of the recombinational signal).

If you are curious: if a third argument is given to biomc2.summarise with a tree file it will compare the trees in this file with the posterior trees and output the corresponding index (tree ID) in the posterior samples (I used this obscure option when comparing with simulations, where we know the true tree and want to see if the sampler can find it). Don't worry bout this option.


### biomc2.distance 
Usage: biomc2.distance tree_file number_of_SPR number_of_replicates > output_file
This program may be run independently from the analysis. Given an initial tree it simulates new topologies whose recombination distance from the previously simulated  tree is given by the user. The initial tree is sampled from 'tree_file', and the program will create a chain of topologies such that all have the same recombination distance from their neighbors. The program may not simulate perfectly the given recombination distance since it applies the specified number of SPRs on the current tree to generate the next tree, but for large numbers we may have cycles. The program tries to minimize this by choosing the prune and regraft edges without replacement. The output_file will contain the tree with uniformly sampled branch lengths, preceded by : the simulation number, real and estimated recombination distances, respectively.
	
The program will also print to the stderr the mean and standard deviation of a few statistics: the overall error, given by the overestimated distance minus the underestimated distance; the modular error, given by the overestimated distance plus the underestimated distance; the overestimation error; and the underestimation error.

### biomc2.sprdist
Usage: biomc2.sprdist tree1 tree2
Calculates d_SPR between a pair of tree files and outputs the histogram of distances.  A file "recomb_leaves.txt" will contain the frequency of recombination for each leaf.  Remember that this frequency is an approximation. (this program is not described in  the paper and we don't use it yet). 



## License
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2008-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)
