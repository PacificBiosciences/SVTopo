use clap::Parser;
use chrono::Datelike;
use std::path::PathBuf;


#[derive(Clone, Parser)]
#[clap(author, version, about, 
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()))]
pub struct Arguments {
    /// pbmm2-aligned BAM filename.
    #[clap(required = true)]
    #[clap(long = "bam")]
    #[clap(value_name = "BAM")]
    pub bam_filename: PathBuf,

    /// Output directory path
    #[clap(required = true)]
    #[clap(long = "svtopo-dir")]
    #[clap(value_name = "STRING")]
    pub outdir: String,

    /// Sample or project ID. No underscores allowed.
    #[clap(required = true)]
    #[clap(long = "prefix")]
    #[clap(value_name = "STRING")]
    pub prefix: String,

    /// (Recommended) structual variant VCF filename.
    /// Requires `--variant-readnames`
    #[clap(required = false)]
    #[clap(long = "vcf")]
    #[clap(value_name = "VCF")]
    pub vcf_filename: Option<PathBuf>,

    /// (Recommended) json with readnames for variant IDs from VCF.
    /// Requires `--vcf`
    #[clap(required = false)]
    #[clap(long = "variant-readnames")]
    #[clap(value_name = "JSON")]
    pub json_filename: Option<PathBuf>,

    /// BED file of regions to exclude from analysis. GZIP files allowed.
    #[clap(required = true)]
    #[clap(long = "exclude-regions")]
    #[clap(value_name = "BED")]
    pub exclude_regions_path: String,

    /// Flag to disable coverage depth and mapping quality filters. Use with caution as many low-confidence images may result.
    #[clap(long = "no-filter", hide = true)]
    pub nofilter: bool,

    /// Flag to output results in unzipped format
    #[clap(long = "write-unzipped", hide = true)]
    pub write_unzipped: bool,

    /// Flag to include connections found using reads that are not haplotagged. This may results in lower-quality outputs.
    #[clap(long = "allow-unphased", hide = true)]
    pub allow_unphased: bool,

    /// Filter threshold for maximum coverage, to remove regions with coverage spikes due to e.g. alignment issues
    #[clap(required = false)]
    #[clap(long = "max-coverage")]
    #[clap(value_name = "INT")]
    #[clap(default_value = "300")]
    pub max_coverage: u32,


    /// Target coordinate in standard {chrom}:{pos}-{end} format
    #[clap(required = false)]
    #[clap(hide = true)]
    #[clap(long = "target-region")]
    pub target_region: Option<String>,

    /// Optional flag to print verbose output for debugging purposes.
    #[clap(long = "verbose")]
    pub verbose: bool,
}

pub fn get_args() -> Arguments {
    Arguments::parse()
}