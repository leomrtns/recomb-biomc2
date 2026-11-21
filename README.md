# biomc2
**biomc2** infers recombination breakpoints in a DNA alignment from changes in phylogenetic topology. The sampler uses a hierarchical Bayesian model with a prior on the distance between neighbouring segment trees. It runs two Metropolis-coupled MCMC chains and reports posterior samples of the mosaic of trees along the alignment.

The model is described in these publications.

- Leonardo de Oliveira Martins, Élcio Leal and Hirohisa Kishino (2008). Phylogenetic Detection of Recombination with a Bayesian Prior on the Distance between Trees. *PLoS ONE* 3(7), e2651. [doi:10.1371/journal.pone.0002651](https://doi.org/10.1371/journal.pone.0002651)
- Leonardo de Oliveira Martins and Hirohisa Kishino (2010). Distribution of distances between topologies and its effect on detection of phylogenetic recombination. *Annals of the Institute of Statistical Mathematics* 62, 145–159. [doi:10.1007/s10463-009-0259-8](https://doi.org/10.1007/s10463-009-0259-8)

Version 2 supports trees spanning multiple 64-bit words and offers RF and Hdist alongside the default approximate dSPR. The likelihood calculation retains the HKY model integrated over exponentially distributed branch lengths.
## Requirements
| Component | Requirements |
| --- | --- |
| Compile the four C executables | A C99-capable compiler, Make and the system maths library |
| Build from a Git checkout or GitHub's automatic source archives | Autoconf 2.70 or later and Automake, in addition to the compiler |
| Run the regression tests | Python 3.10 or later |
| Run `biomc2.plot` | **R**, with `Rscript` on `PATH`. The script uses base R and requires no additional R packages. |

The C executables are `biomc2`, `biomc2.summarise`, `biomc2.distance` and `biomc2.sprdist`. The fifth command, `biomc2.plot`, is an R script installed alongside them. R is needed when drawing figures. The sampler and C summariser run independently of R. If R is available during `make check`, the tests also generate a PDF.

Linux and macOS are the intended platforms. The GitHub workflow builds and tests the Conda package on both. Individual builds should be checked with `make check` before use.
## Compile and install
### From Git
```sh
git clone git@github.com:leomrtns/recomb-biomc2.git
cd recomb-biomc2
autoreconf -fi
mkdir build
cd build
../configure --prefix="$PWD/install"
make -j4
make check
make install
export PATH="$PWD/install/bin:$PATH"
```

Choose an installation prefix suitable for your system. Set `CC`, `CFLAGS`, `CPPFLAGS` or `LDFLAGS` when configuring to select a compiler or supply build flags. For example, `../configure CFLAGS="-O3 -g" --prefix="$PWD/install"` requests an optimised build with debugging information. The default uses the compiler flags chosen by Autoconf, normally `-g -O2`.

`make install DESTDIR="$staging_directory"` supports staged installation. The optional `--enable-static-binary` requests static linking where the toolchain supplies the required libraries. The optional `--enable-optimization` selects `-O3` when the user has supplied no `CFLAGS`.
### From a distribution tarball
A tarball made with `make dist` includes `configure` and the generated Makefiles. After extracting it, follow the instructions above from `mkdir build`. GitHub's automatically generated ZIP and tar archives contain the Git checkout and therefore require `autoreconf -fi`.

Maintainers can run `make distcheck` to check compilation, tests and installation from a distribution tarball.
### Conda
A local recipe is included in `conda/`. In an environment with conda-build installed, run the following command from the repository root.

```sh
conda build --override-channels -c conda-forge conda
```

The recipe takes the current source directory, uses Conda's C compiler and installs the executables under `$PREFIX/bin`. It declares R as a runtime dependency because the package includes `biomc2.plot`. Its tests run the installed programs on the supplied example. The recipe can be built locally or in CI. Publication to Bioconda is a separate step.
## Run the example
The supplied `example/` contains the alignment, initial tree, control file and a historical posterior sample. Work in a fresh directory so that each run has its own output files. From the repository root, with the installed commands on `PATH`, run the following.

```sh
mkdir example-run
cd example-run
cp ../example/paup.nexus ../example/paup.tre ../example/ctrl.biomc2 .
biomc2 ctrl.biomc2
biomc2.summarise post.tree post.dist
biomc2.plot breakpoint-profile.tsv recomb_freq.pdf
```

The last command runs the R script. Check `Rscript --version` if the shell cannot find its interpreter. To call the source script explicitly, use `Rscript ../scripts/biomc2.plot breakpoint-profile.tsv recomb_freq.pdf`.

To draw the figure from the historical posterior sample without running MCMC, replace the sampler and summariser commands with `biomc2.summarise ../example/post.tree ../example/post.dist`. The final plotting command is the same.
## Sampler input
`biomc2 control_file` reads an aligned DNA dataset and one or more initial trees. Filenames in the control file are resolved relative to the directory in which the program runs. The alignment and every tree must have the same set of taxon labels.

The alignment uses NEXUS format with a `MATRIX` block and the `INTERLEAVE` keyword. A single block containing the whole alignment is accepted. The sampler requires at least four taxa and two aligned sites. Gaps and recognised IUPAC ambiguity codes contribute the corresponding sets of possible states.

The tree file uses NEXUS format. Trees must be fully bifurcating, with a trifurcation permitted at the root for an unrooted representation. A rooted representation is converted to an unrooted tree internally. The sampler selects an initial tree at random if several are provided. Branch lengths, when present, determine an initial average rate. Otherwise the initial rate is estimated from the alignment.
### Control file
The required entries are `treefile` and `seqfile`. Use `name = value`, with comments in NEXUS brackets. Parameter names are case-insensitive. Filenames may contain no whitespace and are limited to 511 bytes. An example minimal control file follows.

```text
treefile = paup.tre
seqfile = paup.nexus
metric = dspr
maxdist = 5
ngen = 10000
nsamples = 100
warmup = 100
burnin = 1000
```

| Parameter | Meaning and default |
| --- | --- |
| `metric` | Tree distance, `dspr` by default. Alternatives are `rf` and `hdist`. |
| `maxdist` | Largest permitted distance between neighbouring trees, in units of the selected metric. The default is the number of taxa minus three. Zero requests the historical binary same/different model. |
| `nsegments` | Number of segments at whose boundaries breakpoints may occur. If omitted, multiple NEXUS `CHARSET` blocks determine the segments. Otherwise the alignment is divided into approximately half as many segments as sites, with a minimum of two. |
| `ngen` | Number of sampling iterations after warm-up and burn-in, default 10000. |
| `nsamples` | Requested number of posterior draws, default 100. Integer thinning determines the actual number, which is recorded in the trace header. |
| `warmup` | Iterations per simulated-annealing warm-up cycle, default 100. There are five cycles. |
| `burnin` | Iterations per burn-in cycle, default 1000. There are five cycles. |
| `prior` | Set to 1 to sample the prior, default 0. |
| `rate`, `kappa` | Means of exponential hyperpriors for average rate and transition/transversion ratio, both default 1. |
| `distance` | Shape and rate of the gamma hyperprior on the segment distance parameter, default `1 1`. |
| `penalty` | Shape and rate of the gamma hyperprior on the distance penalty, default `1 1`. |

The supplied control files also document proposal tuning. `update` supplies four perturbation sizes for rate, kappa, distance parameter and penalty. `proposal` controls the interval between SPR proposals and the number of extra topology update cycles. `annealing` sets the final inverse temperature for warm-up. `mc3` specifies the heated-chain inverse temperature and swap interval. `awadhi` and `minisampler` each take an inverse-temperature factor and a number of internal cycles.
### Tree distances
| Value | Calculation | Interpretation |
| --- | --- | --- |
| `dspr` | Approximate unrooted SPR distance, adapted from the splitset implementation in genefam-dist | The original model's approximation to the number of recombination events between segment trees |
| `rf` | Robinson–Foulds symmetric difference, counting splits present in just one tree | A topology discrepancy in split units. One NNI gives RF = 2. |
| `hdist` | Minimum assignment cost between unmatched splits. A split-pair cost is the smaller Hamming distance across the two orientations. | A discrepancy in leaf-membership units |

RF avoids the assignment calculation. Hdist solves one assignment. Approximate dSPR can require several assignment and pruning rounds. The best choice depends on the size and similarity of the trees. Changing the metric changes the prior and can change the posterior distribution.

Set `maxdist` in the selected distance's units. For example, `metric = rf` with `maxdist = 1` excludes all topology changes. Values above the taxon count can be meaningful for Hdist. RF and Hdist quantify topological discrepancy and require that interpretation when reporting results.

Multiword splits remove the original machine-word limit. Partial likelihoods are rescaled to prevent numerical underflow on large trees. Runtime still depends on the number of taxa and distinct site patterns, as well as the MCMC settings. The likelihood and distance implementations retain working memory between proposals.
## Posterior files
The sampler writes `post.tree` and `post.dist` in the working directory. A prior run writes `prior.tree` and `prior.dist`. Use separate working directories for independent runs.

`post.tree` is a NEXUS file whose trees follow the order of draws and segments in the distance trace. `post.dist` starts with the number of segments, number of stored draws and metric name. Its second line gives boundary positions and the alignment length. Every remaining line contains the iteration, log score, total tree distance, number of breakpoints, prior mean rate and prior kappa. Each `|` introduces a boundary index and the distance at that boundary. The boundary index is zero-based and refers to the second header line.

The ancestral names `nSPR` and `dSPR` remain in parts of the C implementation. The values in a trace use the metric recorded in its header. In the binary model the header records `Metric binary`. The summariser reads historical traces without a metric field as dSPR.
## Summarisation and figures
### `biomc2.summarise`
```sh
biomc2.summarise post.tree post.dist
```

The C summariser reads the posterior and estimates a centroid mosaic and a mosaic based on piecewise breakpoint quantiles. It uses the posterior weights of duplicate draws when calculating probabilities. It also compares neighbouring MAP trees with dSPR to identify candidate recombinant leaves.

| Output | Contents |
| --- | --- |
| `mosaicTreeFreq.txt` | Segment boundaries, frequencies and IDs of the two most frequent trees, estimated mosaic boundaries and breakpoint probability |
| `mosaicTrees.tre` | Relevant topologies, with tree IDs used in `mosaicTreeFreq.txt` |
| `mosaicTrees-MAP.tre` | Topologies along the segment-wise MAP mosaic |
| `recSequences.txt` | Candidate pruned leaves at changes in the MAP mosaic and their totals |
| `breakpoint-profile.tsv` | Posterior weights, breakpoint counts and distance sums for plotting or further analysis |

Pruned leaves depend on an approximate SPR reconstruction. Several sequences of moves can explain the same pair of trees.

The centroid is the sampled mosaic with minimum weighted total distance to the other sampled mosaics. Ties are resolved using the sampled total tree distance. The piecewise quantile mosaic divides pooled breakpoint mass into as many equal-mass partitions as the modal number of breakpoints, then selects a median boundary for each partition.

The nine columns of `mosaicTreeFreq.txt` are `loc`, `mapfreq`, `maptre`, `map2freq`, `map2tre`, `mapBP`, `mosBP`, `medBP` and `recfreq`. `loc` is the last site of the segment. Tree IDs are one-based. `mapBP` indicates a change between the MAP tree for this segment and the next. The centroid and piecewise median tree IDs appear at their mosaic boundaries in `mosBP` and `medBP`; zero denotes a boundary absent from that estimate. `recfreq` is the posterior probability of a breakpoint at `loc`.

An optional third argument gives a NEXUS file of guide trees. The summariser prints their indices among the posterior topologies, which is useful when comparing inference with simulated truth.
### `biomc2.plot` — base R
```sh
biomc2.plot [breakpoint-profile.tsv] [recomb_freq.pdf]
```

Both arguments are optional. The R script produces a PDF with two panels. The first uses all draws; the second uses draws with the modal breakpoint count. Grey lines show posterior mean tree distance. Blue points show breakpoint probability. Orange bars summarise the central 95% of each equal-mass partition of pooled breakpoint locations, with a diamond at its median. Grey rug marks identify candidate boundaries with no observed breakpoints in that panel.

These intervals describe the pooled breakpoint signal within each partition. Their interpretation depends on how well the posterior separates individual breakpoint regions.

The TSV contains the columns below and can also be read directly by other analysis software.

| Column | Meaning |
| --- | --- |
| `group` | `all` or `modal` |
| `position` | Alignment boundary, expressed as the last site to its left |
| `weight` | Number of original posterior draws in the group |
| `break_weight` | Number of those draws with a breakpoint at this boundary |
| `distance_sum` | Sum of the sampled distances at this boundary |
| `modal_count` | Modal number of breakpoints |
| `alignment_length` | Number of aligned sites |
| `metric` | Metric recorded in the sampler trace |

Divide `break_weight` by `weight` for breakpoint probability, or `distance_sum` by `weight` for posterior mean distance. Conditional groups use their own total weight. Every candidate boundary is included, even if its breakpoint count is zero.
## Distance utilities
`biomc2.sprdist tree1 tree2 [dspr|rf|hdist]` compares two NEXUS sets of trees on the same labelled tips. It prints a distance histogram weighted by observed tree frequencies. The default is dSPR, which also writes candidate recombinant-leaf frequencies to `recomb_leaves.txt`.

`biomc2.distance tree_file number_of_replicates` generates trees by applying random SPR moves. It reports the number applied, estimated dSPR and half the RF symmetric difference. The resulting trees are written to `spr.trees`.
## Development checks
`make check` exercises the example posterior and all three sampler metrics, including trees across the 64-bit boundary. It compares cached and full likelihoods against an independent reference calculation, including accepted and rejected proposals. With R on `PATH`, it also generates a figure.

For a debug build with GCC or Clang, use `CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"` and the same sanitiser flags in `LDFLAGS` when configuring. Numerical rescaling is required for large-tree likelihoods. Optimisation flags that discard IEEE floating-point semantics, such as `-ffast-math`, need separate validation.
## Licence
The C sources declare GPL version 2 or later, and `COPYING` contains GPL version 2. The adapted split and assignment code retains its original copyright notices.
