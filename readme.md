# HiFi-SVTopo 

SVTopo is a software tool designed to improve visual review of complex structural variation using PacBio HiFi data. SVTopo scans through aligned HiFi reads (whole-genome BAM or partial/BAMlet data) to identify breakends, outputs connected breakends as JSON data, and plots the connected complex structural rearrangements as high-quality images:

![](docs/imgs/complex_fully_connected.png)

HiFi-SVTopo is distributed as two interdependent tools, both necesary for the two-step image generation process:
- HiFi-SVTopo (`svtopo`): a Rust binary for fast BAM parsing into small JSON files. `svtopo` is pronounced "es-vee-tope-oh"
- HiFi-SVTopoVz (`svtopovz`): a Python plotter for the JSON data. `svtopovz` is pronounced "es-vee-tope-oh-viz"

## Documentation
* [Installation](docs/installation.md)
* [HiFi-SVTopo](docs/svtopo_usage.md)
* [HiFi-SVTopoVz usage](docs/svtopovz_usage.md)
* [Results and interpretation](docs/result_interpretation.md)


## Support information
SVTopo is a pre-release software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that SVTopo lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As SVTopo is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any HiPhase release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
