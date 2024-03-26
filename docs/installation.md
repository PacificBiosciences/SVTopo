## Installation 

The `svtopo` Rust utility is available from the [Releases](https://github.com/PacificBiosciences/HiFi-SVTopo/releases) page.

It can be downloaded and run directly on Linux systems, as described in the **Quick Start** section of [HiFi-SVTopo usage](svtopo_usage.md)


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
SVTopoVz v0.3.1
usage: svtopovz [-h] [-v] -j JSON -p OUT_PREFIX [-b BED] [--image-type {png,jpg,jpeg,svg,pdf}] [--log-level {debug,info}]

options:
  -h, --help            show this help message and exit
  -v, --version         Installed version (0.3.1)
  -j JSON, --json JSON  path to json-formatted complex SV data, extracted from BAM file using SVTopo. GZIP allowed. (default: None)
  -p OUT_PREFIX, --out-prefix OUT_PREFIX
                        prefix for output files. May include a directory path. (default: None)
  -b BED, --bed BED     paths to genome annotations in BED file format (default: None)
  --image-type {png,jpg,jpeg,svg,pdf}
                        type of image to generate (default: png)
  --log-level {debug,info}
                        set log level (default: info)
```
