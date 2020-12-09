
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NanoBAC R package

<!-- badges: start -->

[![Travis build
Status](https://travis-ci.org/pgpmartin/NanoBAC.svg?branch=master)](https://travis-ci.org/pgpmartin/NanoBAC)
[![Appveyor build
status](https://ci.appveyor.com/api/projects/status/d55m39h9fk8ijxtc?svg=true)](https://ci.appveyor.com/project/pgpmartin/nanobac)
[![codecov code
coverage](https://codecov.io/gh/pgpmartin/NanoBAC/branch/master/graph/badge.svg)](https://codecov.io/gh/pgpmartin/NanoBAC)
<!-- badges: end -->

The NanoBAC R package provides functions for the assembly of sequences
from BAC or other large plasmids obtained using long read sequencing
technologies such as [Oxford Nanopore](https://nanoporetech.com/) or
[Pacific Biosciences](https://www.pacb.com/).

## Installation

<!-- 
You can install the released version of NanoBAC from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("NanoBAC")
```
-->

You can install the development version from
[GitHub](https://github.com/pgpmartin/NanoBAC) with:

``` r
# install.packages("devtools")
  devtools::install_github("pgpmartin/NanoBAC")
```

A [Singularity](https://sylabs.io/singularity/) container with the
NanoBAC package and its dependencies is also available (see definition
on
[github](https://github.com/pgpmartin/SingIMG/blob/master/Singularity.R36_NanoBAC))  
First build the container `nanobac.sif`using: `singularity build
nanobac.sif shub://pgpmartin/SingIMG:r36_nanobac`

To open a shell inside the container, use: `singularity shell
nanobac.sif`  
To start the R version from the container, use: `singularity exec
nanobac.sif R`  
To run a script that uses the NanoBAC package, use: `singularity exec
nanobac.sif Rscript myRscript.R`

<!-- Add info about singularity-hub image -->

## Assembling repetitive regions of the genomes

Sequencing large plasmids using a long read technology is especially
useful to obtain high quality sequences of *repetitive regions*, a
prerequisite to identify small variations in these regions and assemble
them in order to close the gaps that still remain in sequenced genomes.
In particular, I developed the NanoBAC package while sequencing BACs
containing repeats of the 45S ribosomal DNA (rDNA) gene of the plant
model Arabidopsis thaliana. In the A. thaliana genome, 45S genes are
organized as tandem repeats of ~10kb each and are located in 2 regions
called Nucleolar Organizer Regions (NORs) located at the end of
chromosome 2 and chromosome 4.

While short read sequencing platforms can identify some variations (such
as SNPs or small indels \<\< read length), they can’t reveal how these
variations combine in individual rDNA gene copies or how these copies
are organized along the NORs. Conversely, long reads from genomic DNA
are still too noisy to reliably identify small variations.

Sequencing BACs, each containing several copies of rDNA genes or other
types of repeats, using long read technologies allows to get the best of
both worlds:

  - obtaining **multiple independent reads** for each BAC and aligning
    them allows to correct the random errors inherent to these long
    reads and to gain high confidence in the identification of all
    variations
  - **long reads** spanning several genes/repeats provide information on
    the combination of small variations present in each copy and on the
    organization of these copies relative to one another

By acquiring sequences for a number of long DNA fragments, it may be
possible to assemble the whole NORs and thus close the biggest gap still
present in the A. thaliana genome, and many other genomes as well. This
is also the necessary step to gain a full understanding of how these
rDNA copies are organized at the chromatin and 3D level and how their
expression is regulated.

## Long reads generated during sequencing

A plasmid or BAC is composed of:

  - a **vector (V)** of known sequence  
  - a **DNA insert (D)** of unknown or partially known sequence

The DNA insert is cloned into the vector using a unique
restriction/cloning site in the plasmid.

The different approaches to plasmid/BAC sequencing on long read
platforms can generate different types of reads that are useful to
assemble a high quality sequence:

  - **VDV reads** start by a piece of the **V**ector sequence, followed
    by the full **D**NA insert sequence and end with a fragment of the
    **V**ector sequence
  - **DVD reads** start by a fragment of the **D**NA insert sequence,
    followed by the full **V**ector sequence and end with a piece of the
    **D**NA insert sequence
  - **VD reads** (or **DV reads**) contain a fragment (or the full)
    **V**ector sequence and a fragment (or the full) **D**NA insert
    sequence

We consider that the DNA insert sequence is unknown and repetitive.
Thus, reads only containing the **D**NA insert sequence but no vector
sequence (“**D reads**”) are generally not useful here because we don’t
where to locate them with confidence along the repeats. In other words,
the known and unique vector sequence plays the important role of
anchoring the alignment between the reads so that we can correct the
random sequencing errors and obtain a high quality assembly of the DNA
insert.

## Functions in the NanoBAC package

The functions in the NanoBAC package are used for the following tasks:

1.  Create an interface between the alignment programs Blast/minimap2
    and R:  
    `readBlast`, `blaST2GR`, `SelectSingularBlastALN`, `read_paf`,
    `paf2gr`

2.  Annotate the reads (especially in terms of read type: VDV, DVD,
    etc.):  
    `AnnotateBACreads`, `getVDVnames`, `getDVDnames`, `getVDnames`

3.  Select specific categories or reads / filtering the reads:  
    `FilterBACreads`, `estimateBACsizeFromVDV`, `selectVDVreads`,
    `splitDVDreads`

4.  identify the junctions between the vector sequence and the DNA
    insert sequence in VDV reads:  
    `findVDVjunctions`

5.  Obtain consensus sequences from the files generated by multiple
    sequence alignment softwares:  
    `consensusFromMSF`, `consmat2seq`

In addition, the `makeVarseq` function allows to introduce random
variations in a DNA sequence in order to create simulated reads.

These function are at the core of the [NanoBAC
pipeline](https://github.com/pgpmartin/NanoBAC_pipeline). Below are
examples of how to use the functions of the NanoBAC R package.  
Start by loading the package:

``` r
library(NanoBAC)
```

## Importing alignment results

In order to characterize the reads obtained during a sequencing run, it
is necessary to align the vector sequence and potential known sequences
that are expected in the DNA insert (e.g. the 18S rDNA sequence). This
can be done efficiently using
[Blast](https://www.ncbi.nlm.nih.gov/books/NBK279690) which, in command
line, has an option to output its result in a convenient tabular form:
`--outfmt=6`.

We can import such data using:
