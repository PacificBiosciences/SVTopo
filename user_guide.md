## Overview
SVTopo creates images to represent chimeric/split alignment evidence for structural variation, with an emphasis on complex variants. For the purposes of SVTopo, these are defined as SVs large enough to for the PBMM2 aligner to align reads in in multiple chimeric alignments. Novel insertions are not supported and by default simple deletions and duplications are omitted from results. Complex SVs in SVTopo therefore consist of inversions, translocations, and multi-breakend combinations of simple rearrangements.


__Processing steps:__
1. Extract reads flagged with the BAM split alignment tag (`SA`)
2. Identify coordinates where alignment clipping sites from multiple reads are clustered
3. Connect pairs of coordinates that are joined by alignments
4. Further connect these pairs to generate graphs of breakends that define structural rearrangements
5. Emit graph representations as a JSON file


## Result interpretation

## Getting started 

### HiFi-SVTopo
The `svtopo` Rust utility is available from the [Releases](https://github.com/PacificBiosciences/HiFi-SVTopo/releases) page.

It can be downloaded and run directly on Linux systems, as described in the **Quick Start** section of [HiFi-SVTopo usage](svtopo_usage.md)


### HiFi-SVTopoVz
The `svtopovz` utility can be downloaded by cloning this reposity with git or by downloading the `Source code (zip)` or `Source code (tar.gz)` assets from the [Releases](https://github.com/PacificBiosciences/HiFi-SVTopo/releases) page.

It can then be installed from source as shown below. It is recommended that this tool be installed in a fresh [Conda](https://conda.io/projects/conda/en/latest/index.html) environment with Python=3.10. 

Installation from source:
```bash
git clone https://github.com/PacificBiosciences/HiFi-SVTopo.git

cd HiFi-SVTopo/HiFi-SVTopoVz/
python setup.py install
```

Test that the installation succeeded:
```bash
svtopovz -h
```
If successfully installed, this will print out the tool version and command-line options: 
```text
svtopovz -h

SVTopoVz v0.1.0
usage: svtopovz [-h] [-v] --json JSON --out-prefix OUT_PREFIX [--bed BED] [--image-type {png,jpg,jpeg,svg,pdf}] [--max-gap-size-mb MAX_GAP_SIZE_MB] [--verbose] [--include-simple-dels]
                [--include-simple-dups]

options:
  -h, --help            show this help message and exit
  -v, --version         Installed version (0.1.0)
  --json JSON           path to json-formatted complex SV data, extracted from BAM file using SVTopo. GZIP allowed. (default: None)
  --out-prefix OUT_PREFIX
                        prefix for output files. May include a directory path. (default: None)
  --bed BED             paths to genome annotations in BED file format (default: None)
  --image-type {png,jpg,jpeg,svg,pdf}
                        type of image to generate (default: png)
  --max-gap-size-mb MAX_GAP_SIZE_MB
                        maximum gap size to show in one panel, in megabases. Default is 0.5mB (default: 0.5)
  --verbose             print verbose output for debugging purposes (default: False)
  --include-simple-dels
                        does not skip simple deletions (defined as having exactly two forward spanned blocks, one forward unspanned) (default: False)
  --include-simple-dups
                        does not skip simple duplications (defined as having exactly two forward spanned blocks, one reverse unspanned) (default: False)

```
## Outputs

## Usage details
