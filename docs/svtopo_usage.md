## SVTopo BAM parsing

### Getting started
**Input:** SVTopo requires phased HiFi BAM input. Phasing with [HiPhase](https://github.com/PacificBiosciences/HiPhase) is strongly recommended and results may be negatively impacted if a different phasing tool is used.

```
Usage: svtopo [OPTIONS] --bam <BAM> --exclude-regions <BED>

Options:
      --bam <BAM>                 pbmm2-aligned BAM filename
      --vcf <VCF>                 (Recommended) structual variant VCF filename. Requires `--variant-readnames`
      --variant-readnames <JSON>  (Recommended) json with readnames for variant IDs from VCF. Requires `--vcf`
      --exclude-regions <BED>     BED file of regions to exclude from analysis. GZIP files and file URLs allowed
      --max-coverage <INT>        Filter threshold for maximum coverage, to remove regions with coverage spikes due to e.g. alignment issues [default: 300]
  -h, --help                      Print help
  -V, --version                   Print version

Copyright (C) 2004-2024     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.
```

## Recommended usage
SVTopo gets the best results when paired with the PacBio [Sawfish](https://github.com/PacificBiosciences/sawfish) SV caller. Sawfish is run first with the [`--report-supoorting reads` option](https://github.com/PacificBiosciences/sawfish/blob/main/docs/user_guide.md#phasing), which generates a JSON file of variant IDs keyed to the readnames for supporting alignments (`supporting_reads.json`). SVTopo accepts both the SV VCF and supporting_reads.json file as inputs, which allow SVTopo to use sawfish SV call coordinates for SV block identification. 

### Example SVTopo workflow with sawfish variant calls
```bash
svtopo\
    --bam HG002_hg38.bam\
    --json-out svtopo_hg002_output.json\
    --variant-readnames sawfish_HG002_supporting_reads.json\
    --vcf HG002_sawfish.vcf.gz \
    --exclude-regions https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz
```

### Example SVTopo workflow without sawfish variant calls
SVTopo can also attempt to identify variants _de novo_ (i.e. without sawfish calls). SVTopo is not intended to be a fully supported variant caller and thus results will be less accurate than with help from sawfish, but this provides a simpler interface:
```bash
svtopo\
    --bam HG002_hg38.bam\
    --json-out svtopo_hg002_output.json\
    --exclude-regions https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz
```

## Runtime and output
SVTopo was benchmarked with a 30x HiFi genome, with and without the recommended sawfish-based option:

* Sawfish benchmark:
      * Runtime 13mins:33secs
      * 3.3 GB RAM
      * 202 complex SV images
* Non-sawfish benchmark:
      * Runtime 10mins:2secs
      * 3.4 GB RAM
      * 216 complex SV images

## Algorithm notes
* Clipped alignments: SVTopo uses chimeric/split alignments to identify signals of structural variation. These are defined as alignments with at least 100 bases of soft-clipping on either end of the alignment. Alignments with MAPQ < 20 are omitted.
* Finding break clusters with alignments: Locations of genomic breaks are identified using alignment clipping locations that are clustered together. These must be within a 10 bp confidence interval of each other, allowing for small differences of alignment and minimal sequencing error. Break cluster positions are given from the middle of the clip sites among alignments within 10 bp.
* Finding break clusters with sawfish output: If a VCF and `supporting_reads.json` are provided from [sawfish](https://github.com/PacificBiosciences/sawfish), these are used to identify additional break locations (the VCF SV POS and END locations). In this case, the alignments from a read assigned to a given variant ID in the sawfish output are compared to the variant breakend locations and the clipping sites nearest the variant break locations are used to assign specific split alignments from the read to the coordinate of the breakend.
* Phase filter: break clusters are only included if the phase is unambiguous. That is defined as having exactly one haplotype ID and one phaseset ID, with at least 2 phased alignments. Unphased alignments may also be included.
* Break connections: Once breaks are identified, they can be connected in pairwise fashion by alignments that are shared between them, or by using VCF entry connections (via VCF record POS/END). Alignment-based connections must have a minimum of two such shared alignments and break connections must be within 1 mb of each other (if on the same chromosome).
* Phased connection of clusters: In some cases phasing is used to connect break locations instead of directly using spanning alignments. In these cases, the clusters must be within 500kb of each other, have the same phaseset ID, and be on the same haplotype.
* Ambiguous sample order: In some cases it may be impossible to determine the order of some genomic blocks. These are given the sample sample_order_index entry in JSON and plotted in images with red outlines.
* The following filters are applied:
  * Coverage of 300x or less (can be changed using the `max-coverage` option)
  * No more than 5% of reads below MAPQ of 5
