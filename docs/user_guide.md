# SVTopo User Guide

## Table of Contents
* [Overview](#overview)
* [Getting Started](#getting-started)
* [Outputs](#outputs)
* [Usage Details](#usage-details)

## Overview

SVTopo creates images to represent chimeric/split alignment evidence for structural variation, with an emphasis on complex variants. For the purposes of SVTopo, these are defined as SVs large enough to for the [pbmm2](https://github.com/PacificBiosciences/pbmm2) aligner to align reads in in multiple chimeric alignments. Novel insertions are not supported and by default simple deletions and duplications are omitted from results. Complex SVs in SVTopo therefore consist of inversions, translocations, and multi-breakend combinations of simple rearrangements.

Processing steps:
1. Extract reads flagged with the BAM split alignment tag (`SA`)
2. Identify coordinates where alignment clipping sites from multiple reads are clustered
3. Connect pairs of coordinates that are joined by alignments
4. Further connect these pairs to generate graphs of breakends that define structural rearrangements
5. Emit graph representations as a JSON file
6. Generate images/HTML viewer from the JSON file with the SVTopoVz utility

## Getting started 

### Installation

#### Install from Conda
SVTopo and the plotting utility SVTopoVz are both available from [Bioconda](https://bioconda.github.io/) on Linux. Assuming you have already installed [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) and [mamba](https://mamba.readthedocs.io/en/latest/), the following example code creates a new environment with python v3.10, activates that environment, and installs both svtopo and svtopovz (as a dependency) into that environment. Note that these mamba commands work the same with conda, but mamba is recommended due to improved execution speed.

```bash
conda config --add channels bioconda
mamba create -n svtopo "python=3.10"
mamba activate svtopo
mamba install -y svtopo
```

SVTopoVz can also be installed directly on MacOs via conda:
```bash
conda config --add channels bioconda
mamba create -n svtopo "python=3.10"
mamba activate svtopo
mamba install -y svtopovz
```

#### Install from GitHub
##### SVTopo installation
The `svtopo` Rust utility is available from the [Releases](https://github.com/PacificBiosciences/SVTopo/releases) page.

It can be downloaded, extracted, and run directly on Linux systems. For example with version 0.2.0:
```bash
wget https://github.com/PacificBiosciences/SVTopo/releases/download/v0.1.1/svvtopo_v0.2.0-x86_64-unknown-linux-gnu.tar.gz
tar -zxvf svtopo_v0.1.1-x86_64-unknown-linux-gnu.tar.gz
svtopo --version
```
This should will the `svtopo` binary to the local directory and print out the version number.

The `svtopo` tool can also been installed directly from source by downloading and building the source code directly:
```bash
git clone https://github.com/PacificBiosciences/SVTopo.git
cd SVTopo/
cargo build --release
```
The `svtopo` binary will be created at `SVTopo/target/releases/svtopo`.

##### SVTopoVz installation
The `svtopovz` utility can be downloaded by cloning this reposity with git or by downloading a Source code asset from the [Releases](https://github.com/PacificBiosciences/SVTopo/releases) page.

It can then be installed from source as shown below. It is recommended that this tool be installed in a fresh [Conda](https://conda.io/projects/conda/en/latest/index.html) environment with Python=3.10. 

Install from source, then test that the installation succeeded:
```bash
git clone https://github.com/PacificBiosciences/SVTopo.git
cd SVTopo/SVTopoVz/
conda create -n svtopo "python=3.10"
conda activate svtopo
python setup.py install
svtopovz -h
```
If successfully installed, this will print out the tool version and command-line options.

### Run tests
For a more complete test of your installation, run the test cases included in the [test/](https://github.com/PacificBiosciences/SVTopo/tree/main/test) directory.
The test script requires the absolute path to your downloaded `svtopo` binary. The following example assumes the binary is located in the `$HOME/bin/` directory:
```bash
cd SVTopo/
bash test/scripts/run_end_to_end_tests.sh $HOME/bin/svtopo_x86_64
```

### Analysis Steps
SVTopo finds and plots SV images in two steps:
1. `svtopo` - Parses the phased BAM, exclude regions BED, and optional sawfish variant call files (VCF and `supporting_reads.json`) to create a complex variant graph representation in JSON format.
2. `svtopovz` - Reads the svtopo output JSON file and uses the complex variant graph information to create plots of the complex variants.

### Example workflow
To plot complex SVs, run SVTopo on the mapped and phased sample BAM, and then run SVTopoVz on the resulting JSON file.

The following example shows how this is done for a sample BAM called `HG002_hg38.bam` with sawfish VCF and `supporting_reads.json` inputs. The sawfish inputs may be omitted, but their inclusion improves results.

```bash
wget https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz
svtopo\
    --bam HG002_hg38.bam \
    --json-out svtopo_hg002_output.json.gz \
    --variant-readnames sawfish_HG002_supporting_reads.json \
    --vcf HG002_sawfish.vcf.gz \
    --exclude-regions cnv.excluded_regions.hg38.bed.gz
mkdir images/
svtopovz \
  --json svtopo_hg002_output.json.gz \
  --out-prefix images/hg002_svtopo
```

## Outputs
### SVTopo step
The primary output of the SVTopo step is the JSON file detailing the complex SVs found. This may be viewed by the user directly but is not optimized for human readability.

### SVTopoVz step
The SVTopoVz step generates images for each of the entries in the complex SV JSON file that:
* is not a simple deletion or tandem duplication
* can be fully resolved by SVTopoVz into a complex SV graph

For detailed notes on interpretation of these plots, see [result interpretation](https://github.com/PacificBiosciences/SVTopo/blob/main/docs/result_interpretation.md).

## Usage details
### Expected compute requirements
SVTopo benchmarks with a 30x HiFi genome, with and without the recommended sawfish-based option:

__Sawfish benchmark:__
* SVTopo
  * Runtime 33 min, 54.33 sec
  * 2.3 GB RAM
  * 62% CPU
* SVTopoVz
  * 1 min, 22.62 sec
  * 0.18 GB RAM
  * 78 complex SV images

__Non-sawfish benchmark:__
* SVTopo
  * Runtime 31 min, 59.88 sec
  * 4.7 GB RAM
  * 60% CPU
* SVTopoVz
  * 0 min, 56.11 sec
  * 0.18 GB RAM
  * 73 complex SV images
  
### Algorithm notes
* Clipped alignments: SVTopo uses chimeric/split alignments to identify signals of structural variation. These are defined as alignments with at least 100 bases of soft-clipping on either end of the alignment. Alignments with MAPQ < 20 are omitted.
* Finding break clusters:
    * With alignments: Locations of genomic breaks are identified using alignment clipping locations that are clustered together. These must be within a 10 bp confidence interval of each other (allowing for small differences of alignment). A minimum of two alignments is required to support a cluster as a potential valid break location.
    * With sawfish output: If a VCF and `supporting_reads.json` are provided from [sawfish](https://github.com/PacificBiosciences/sawfish), these are used to identify additional break locations (by using the VCF SV POS and END locations. The specific alignment from a chimeric read is assigned to a coordinate break location by identifying the closest pair of read clipping coordinate and variant breakend coordinate.
* Break connections: Once breaks are identified, they can be connected in pairwise fashion by alignments that are shared between them. They may also be connected by using VCF entry connections (via VCF record POS/END). Alignment-based connections must have a minimum of two such shared alignments and must be within 1 mb of each other (if on the same chromosome).
    * Phased connection of clusters: If a breakend lacks direct connections to another breakend via alignments, SVTOpo searched for breakends up or down-stream for 500kb and connects them if the reads supporting both breaks have the same phaseset ID and are on the same haplotype.
* Ambiguous sample order: In some cases it may be impossible to determine the order of some genomic blocks. These are given the sample sample_order_index entry in JSON and plotted in images with red outlines instead of the standard dark grey.
* The following filters are applied:
  * Coverage of 300x or less (can be changed using the `max-coverage` option)
  * No more than 5% of reads below MAPQ of 5