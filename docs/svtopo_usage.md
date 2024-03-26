## Cnidaria BAM parsing

### Quick Start
**Input:** The `cnidaria` utility requires phased HiFi BAM input. Phasing with [HiPhase](https://github.com/PacificBiosciences/HiPhase) is strongly recommended and results may be negatively impacted if a different phasing tool is used.

```
Usage: cnidaria [OPTIONS] --bam <BAM>

Options:
  -b, --bam <BAM>             pbmm2-aligned BAM filename
  -p, --output-prefix <PATH>  Output prefix, which may begin with a directory path [default: ""]
      --keep-dels             Optional flag to include simple deletions. These are filtered out by default
      --verbose               Optional flag to print verbose output for debugging purposes
  -h, --help                  Print help
  -V, --version               Print version
```

**BAM parsing example**
Execution example (replace paths with your own):
```bash
cnidaria --bam /path/to/reads.bam > /path/to/results.json
```


## Advanced usage
#### Including deletions
By default, `cnidaria` will omit events predicted as being simple deletions. Add the `--keep-dels` flag to retain these records. 

Deletions are defined here as connected events where two spanned blocks occur and are joined by a single unspanned block, all in forward orientation.

#### Setting a JSON prefix
The `--output-prefix` option allows users to pass in an explicit prefix for the JSON output, which may consist of a relative or absolute path and an optional file name prefix for the JSON output. If used, JSON output is written in gzipped format. For example, the following input `--prefix mydirectory/sample_id` will generate a JSON file `mydirectory/sample_id_cnidaria.json.gz`. Note that the `cnidaria_plotter` step accepts either zipped or unzipped JSON inputs.


## Runtime and output
`Cnidaria` was tested with a 30x HiFi genome, in default mode and with the `--keep-dels` option.

### Default mode
* 28m20s run time
* 100 images generated

### `--keep-dels` mode
* 29m20s run time
* 385 images generated