## Snakemake workflow
**This section is under construction** 

Cnidaria requires HiFi BAM data that has been phased with [HiPhase](https://github.com/PacificBiosciences/HiPhase). A Snakemake workflow is therefore provided for enhanced ease-of-use. This workflow performs the phasing preparation step.

### Caveat
This workflow will not support all use cases and the official [PacBio WGS Variant Pipeline](https://github.com/PacificBiosciences/HiFi-human-WGS-WDL) is recommended as a more fully developed workflow, with the `merged_haplotagged_bam` output from that pipeline as input to `Cnidaria`.

### Cnidaria snakemake inputs and environment
This snakemake workflow (in `workflow/snakefile`) requires the following inputs:
* Aligned BAM(s)
* Reference genome (in FASTA format)
* Reference index file (in FAI format)
* The path to a `Cnidaria` Rust binary (from Releases)

Additionally, this must be executed in a Linux environment with `Singularity` and `Conda` installed.

### Cnidaria snakemake execution
The `Cnidaria` snakemake workflow expects all inputs to be stored in a config file, such as the example in `workflow/config.json` file. Replace example values with your own data.

To run the workflow, use the following command, filling in the --cores and --singularity-args with your own correct information:
```bash
snakemake \
    --cores {} \
    --use-singularity \
    --use-conda\
    --singularity-args "--bind {} "\
    --configfile config.json
```

Be aware that proper file system binding with the singularity args option is critical to proper function. Refer to [the Singularity](https://docs.sylabs.io/guides/3.0/user-guide/bind_paths_and_mounts.html) and [Snakemake](https://snakemake.readthedocs.io/en/v5.7.3/snakefiles/deployment.html) documentation for more details.

### Config file
The JSON configuration file is expected to follow the following example structure:
```json
{
    "bam_files": {
        "sample": "/path/myfile.bam"
    },
    "reference": "/path/myref.fasta",
    "reference_idx": "/path/myref.fasta.fai",
    "cnidaria_binary": "/path/cnidaria"
}
```
The following fields are required:
* `bam_files`: a dictionary containing unique sample identifiers mapped to HiFi BAM file paths as shown here by `"sample": "/path/myfile.bam"`.
* `reference`: the path to a reference genome in FASTA format.
* `reference_idx`: the path to a reference index file in FAI format.
* `cnidaria_binary`: the path to a local copy of the `cnidaria` binary.

An example configuration file is included at `workflow/config.json` and may be used as a template.


### Cnidaria workflow overview
The major steps taken in processing input BAMs are as follows:
1. Remove any existing phase tags from the BAM (`PS`/`HP` tags) and index them with [samtools view](http://www.htslib.org/doc/samtools-view.html)
2. Call small variants with [DeepVariant](https://github.com/google/deepvariant)
3. Phase and haplotag the BAMs with [HiPhase](https://github.com/PacificBiosciences/HiPhase)
4. Extract complex SV evidence and generate plots of the results using Cnidaria