## SVTopo BAM parsing

### Getting started
**Input:** The SVTopo utility requires phased HiFi BAM input. Phasing with [HiPhase](https://github.com/PacificBiosciences/HiPhase) is strongly recommended and results may be negatively impacted if a different phasing tool is used.

```
Usage: svtopo [OPTIONS] --bam <BAM> --genome <genome>

Options:
      --bam <BAM>                 pbmm2-aligned BAM filename
      --vcf <VCF>                 (Recommended) structual variant VCF filename. Requires `--variant-readnames`
      --variant-readnames <JSON>  (Recommended) json with readnames for variant IDs from VCF. Requires `--vcf`
      --genome <genome>           Genome selection [possible values: hg38, hg19, other]
      --max-coverage <INT>        Filter threshold for maximum coverage, to remove regions with coverage spikes due to e.g. alignment issues [default: 300]
  -h, --help                      Print help
  -V, --version                   Print version

Copyright (C) 2004-2024     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.
```

**BAM parsing example**
Execution example (replace paths with your own):
```bash
svtopo --bam /path/to/reads.bam > /path/to/results.json
```

## Recommended usage
SVTopo gets the best results when paired with the PacBio [Sawfish](https://github.com/PacificBiosciences/sawfish) SV caller. Sawfish is run first with the `--report-supoorting reads` option, which generates a JSON file of variant IDs keyed to the readnames for supporting alignments (`supporting_reads.json`). SVTopo accepts both the SV VCF and supporting_reads.json file as inputs, which allow SVTopo to use Sawfish SV call coordinates for SV block identification. 

### Example SVTopo workflow
```bash

```

## Runtime and output
SVTopo was benchmarked with a 30x HiFi genome, using the recommended VCF option:

#### Default mode
* 91m42s run time
* 586 SV entries in JSON

#### No-filter mode
* 10m43s run time
* 838 SV entries in JSON

## Algorithm notes
* Clipped alignments: SVTopo uses chimeric/split alignments to identify signals of structural variation. These are defined as alignments with at least 100 bases of soft-clipping on either end of the alignment. Alignments with MAPQ < 20 are omitted.
* Break clusters: Locations of genomic breaks are identified using alignment clipping locations that are the same or nearly the same. These are defined as being within 10 bp of each other, allowing for small differences of alignment and minimal sequencing error. Break cluster positions are given from the middle of the clip sites among alignments within 10 bp.
  * Break clusters may be detected via [Sawfish]() VCF entries using the '--vcf' and `--suporting-reads` command-line options
* Phase filter: break clusters are only included if the phase is unambiguous. That is defined as having exactly one haplotype ID and one phaseset ID, with at least 2 phased alignments. Unphased alignments may also be included.
* Break connections: Once breaks are identified, they can be connected in pairwise fashion by alignments that are shared between them, or by using VCF entry connections (via VCF record POS/END). Alignment-based connections must have a minimum of two such shared alignments and break connections must be within 1 mb of each other (if on the same chromosome).
* Phased connection of clusters: In some cases phasing is used to connect break locations instead of directly using spanning alignments. In these cases, the clusters must be within 500kb of each other, have the same phaseset ID, and be on the same haplotype.
* Ambiguous sample order: In some cases it may be impossible to determine the order of some genomic blocks. These are given the sample sample_order_index entry in JSON and plotted in images with red outlines.
* The following filters are applied:
  * Coverage of 300x or less (can be changed using the `max-coverage` option)
  * No more than 5% of reads below MAPQ of 5
