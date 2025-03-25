[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/svtopo/README.html)

<h1 align="center"><img width="250px" src="docs/imgs/svtopo.svg"/></h1>

<h1 align="center">SVTopo</h1>

<h4 align="center">Complex structural variant visualization for HiFi sequencing data</h3>

SVTopo represents complex structural variants (SVs) to enhance PacBio HiFi data. SVTopo determines connections between SV breakends using aligned reads and, optionally, variant calling output from the [sawfish](https://github.com/PacificBiosciences/sawfish) variant caller. SVTopo outputs connected breakends as JSON data and uses the SVTopoVz python utility to plot complex structural rearrangements as high-quality images.

Key Features (see [Result Interpretation](docs/result_interpretation.md) for more details):
* Top section: illustration of connected aligned sequences as blocks, showing order and orientation in the sample genome relative to the reference
  * Alignment support strength shown by block weight
  * Block sizes encoded into legend
* Optional gene annotations
* Bottom section: chained plot highlighting deletion, duplication, and inversion effects on rearranged genome structure

Example SVTopo image from HG002 rearrangement
![](docs/imgs/complex_fully_connected.png)

## Getting started
* See the [Getting started](docs/user_guide.md#getting-started) section in the [User Guide](docs/user_guide.md) to start using SVTopo
* See also [Result interpretation](docs/result_interpretation.md) for examples and explanation of SVTopo images


## Support information
SVTopo is a pre-release software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that SVTopo lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As SVTopo is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any HiPhase release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.