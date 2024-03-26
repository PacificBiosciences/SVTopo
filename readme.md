# Cnidaria 

Cnidaria is a software tool designed to improve visual review of complex structural variation using PacBio HiFi data. Cnidaria scans through aligned HiFi reads (whole-genome BAM or partial/BAMlet data) to identify breakends, outputs connected breakends as JSON data, and plots the connected complex structural rearrangements as high-quality images:

![](docs/imgs/dup_example.png)

Cnidaria is distributed as two complimentary tools, both necesary for the two-step image generation process:
- `Cnidaria`: a Rust binary for fast BAM parsing into small JSON files
- `Cnidaria_plotter`: a Python plotter for the JSON data

## Documentation
* [Installation](docs/installation.md)
* [Cnidaria usage](docs/cnidaria_usage.md)
* [Cnidaria plotter usage](docs/cnidaria_plotter_usage.md)
* [Results and interpretation](docs/result_interpretation.md)
* [Snakemake workflow](docs/snakemake_workflow.md)


## Support information
Cnidaria is a pre-release software intended for research use only and not for use in diagnostic procedures. 
While efforts have been made to ensure that Cnidaria lives up to the quality that PacBio strives for, we make no warranty regarding this software.

As Cnidaria is not covered by any service level agreement or the like, please do not contact a PacBio Field Applications Scientists or PacBio Customer Service for assistance with any HiPhase release. 
Please report all issues through GitHub instead. 
We make no warranty that any such issue will be addressed, to any extent or within any time frame.

### DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
