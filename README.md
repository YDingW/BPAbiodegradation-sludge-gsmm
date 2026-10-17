# BPA biotransformation in activated sludge: modelling code

Code for the genome-scale metabolic modelling, the MICOM community simulations and the Mfuzz clustering of
"Community-level bisphenol A biotransformation in an activated sludge microcosm". The raw sequencing data are
in NCBI under BioProject [PRJNA1412122](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1412122).

    GSMM/     metabolic models and community simulation, one script per step, and the medium
    mfuzz/    clustering of the MAG time series

Each script is self-contained: its parameters are set at its top, and its header states what it does, the
environment it needs and its input and output files. `GSMM/environment.yml` has the Python environment of
`GSMM/`; the R packages of `mfuzz/` are listed at the top of those scripts.

## GSMM/ — in the order the scripts are run

**growth_rates.py** — growth rate of every MAG at every simulated time point. Over an interval of h hours the
whole community grows at

    r = mean over the three replicates of ln(N_end / N_start) / h        (N = 16S rRNA gene copies)

and the relative abundance a_i of MAG i (metaWRAP) changes with the difference between its own growth rate and
r, so

    mu_i = ln(a_i,end / a_i,start) / h + r

A negative mu_i is set to 0. The time point 0 h uses the interval 0–24 h, the time points 24, 26 and 30 h use
24–30 h.

**eflux_bounds.py** — community members and E-Flux bounds, from the metatranscriptome. The expression of a gene
at a time point is the square root of the mean of its three replicates. A MAG is a member when transcripts are
found for at least 5 % of its protein-coding genes. With G the 99th percentile of the expression of all
detected protein-coding genes of all MAGs over the four time points, the coefficient of a reaction is

    c = min(1, score / G)

where score is its GPR evaluated on gene expression (maximum over OR, minimum over AND). c is not allowed below
the floor of its MAG, the 10th percentile of min(1, expression / G) over the detected protein-coding genes of
that MAG; reactions with a GPR but no data and metabolic reactions without a GPR get the floor. Every non-zero
bound is then multiplied by c, so a reaction keeps its direction:

    lower bound -> c x lower bound (if negative),  upper bound -> c x upper bound (if positive)

The biomass reaction, spontaneous reactions and exchange and transport reactions without a GPR keep their
bounds.

**run_community.py** — the community simulation of one time point. The models of the members, each with its
E-Flux bounds, are merged into a MICOM community with the relative abundances of the members, and the medium of
the time point is applied. Then

1. every member is required to grow at the rate from `growth_rates.py`: maximise t subject to mu_i >= t x
   target_i for every member; if t < 1, the members whose constraint is binding are pinned at the rate they
   reached and the problem is solved again for the others;
2. with the members held at those rates, the fluxes are chosen by expression-weighted parsimonious FBA,

       minimise  sum_j price_j |v_j|,   price_j = (1 / c_j)^2 for a reaction with an expression coefficient,
                                        price_j = 1600 for every other reaction

so that a reaction of a highly expressed gene is cheap and a reaction of a weakly expressed gene is expensive.
The solution is checked: no member below its rate, the carbon leaving the community does not exceed the carbon
entering it, and no flux through boundary reactions outside the medium.

**exchange_tables.py** — the metabolite exchange between the members, from the simulations of the four time
points. For a member and a metabolite, the flux of the member's exchange reaction is the specific flux
(negative = uptake, positive = secretion), and

    weighted flux = specific flux x abundance share of the member     (mmol per gDW of community per h)
    carbon flux   = weighted flux x carbon atoms of the metabolite
    cross-fed     = min(total secretion, total uptake) of the metabolite

A metabolite secreted by several members and taken up by several others goes through a shared pool, so a
donor–receiver pair gets the cross-fed amount in proportion to the donor's share of the secretion and the
receiver's share of the uptake.

**energy_cycles.py** — energy-generating cycles (EGCs) of every model, tested with the bounds the model
actually has. Every boundary reaction is closed, a dissipation reaction for an energy carrier is added (for
example atp + h2o -> adp + pi + h) and its flux is maximised; a flux > 0 means the model can charge that
carrier from nothing. The test follows Fritzemeier et al. 2017, as in MEMOTE.

**memote_summary.py** — summary of the MEMOTE reports of the draft and the curated models: one table with the
metrics of every model, one with the median, mean, minimum and maximum per set, and the distribution of every
metric over the 90 models of each set.

**medium.tsv** — the medium of the community simulations: one row per exchange reaction of the community, one
column per time point, with the maximal uptake of the community in mmol per gDW per h. It is a mineral medium
with BPA (0.0365, the measured biotransformation rate) and oxygen; glucose stands for the other carbon sources
of the sludge and is supplied at 0 and 30 h only. A component with 0 is not part of the medium of that time
point.

## mfuzz/ — clustering of the MAG time series

**mfuzz_clustering.R** — soft clustering of the temporal profiles. Every profile is transformed with
log10(x + 1) and standardised to mean 0 and standard deviation 1 over its sampling points, then clustered by
fuzzy c-means (`Mfuzz::mfuzz`) with c = 3 and the fuzzifier m estimated with `Mfuzz::mestimate`; the random
start is fixed with `set.seed(2026)`. A MAG is a core member of its cluster when its membership is at least
0.5. Run once for the MAG abundances and once for their total transcript abundance.

**trend_plot.R** — trend plot of the clusters: one panel per cluster, one line per MAG coloured by its
membership, the black line the cluster centroid.

**sankey_plot.py** — Sankey diagram of the MAGs across abundance cluster, family and expression cluster;
ribbon width is the number of MAGs.
