use std::collections::HashSet;
use std::time::SystemTime;
use std::{collections::HashMap, path::PathBuf};

use log::{debug, error, info, LevelFilter};
use std::env;
use svtopo::cli::{get_args, Arguments};
use svtopo::cluster_connector::connect_clusters;
use svtopo::containers::{ComplexSVCalls, Connection, Coordinate, EventGraph, TargetCoordinate};
use svtopo::containers::{ExcludeRegions, FwdStrandSplitReadSegment};
use svtopo::event_graph_builder::build_event_graphs;
use svtopo::graph_annotator::annotate_graphs;
use svtopo::ingester::{
    get_read_info_from_json, get_sample_from_bam, get_split_alignments,
    get_split_alignments_from_region, get_vcf_breaks, load_exclude_regions,
};
use svtopo::utils::is_local_file;
use svtopo::{block_filter, cluster_finder, cluster_support_builder, result_writer};

fn set_up() -> (Arguments, Option<TargetCoordinate>) {
    let args = get_args();
    let filter_level: LevelFilter = match args.verbose {
        false => LevelFilter::Info,
        true => LevelFilter::Debug,
    };
    env_logger::builder()
        .format_timestamp_millis()
        .filter_level(filter_level)
        .init();

    let version = env!("CARGO_PKG_VERSION");
    info!("\nRunning HiFi-SVTopo v{version}\n");

    let cmd: Vec<String> = env::args().collect();
    let cmd_str = cmd.join(" ");
    debug!("Run command: {cmd_str}");
    debug!("v{version}\n");

    let mut has_vcf = false;
    let mut has_json = false;
    if args.vcf_filename.is_some() {
        has_vcf = true;
    }
    if args.json_filename.is_some() {
        has_json = true;
    }
    if has_json && !has_vcf {
        error!("`--variant-readnames` json input requires `--vcf` input");
        std::process::exit(exitcode::CONFIG);
    }

    if let Some(region_str) = args.target_region.clone() {
        let target_region = TargetCoordinate::new(region_str);
        return (args, Some(target_region));
    }

    if !is_local_file(&args.exclude_regions_path) {
        error!(
            "Exclude regions file {} not found",
            args.exclude_regions_path
        );
        std::process::exit(exitcode::CONFIG);
    }

    let path = std::path::Path::new(&args.outdir);
    if !path.exists() || !path.is_dir() {
        error!("outdir {} does not exist", args.outdir,);
        std::process::exit(exitcode::CONFIG);
    }
    if args.prefix.contains('_') {
        error!("Prefix does not allow underscores");
        std::process::exit(exitcode::CONFIG);
    }

    (args, None)
}

fn log_time(start_time: SystemTime) {
    let elapsed_time = start_time.elapsed().unwrap().as_secs();
    let hours = elapsed_time / 3600;
    let minutes = (elapsed_time % 3600) / 60;
    let seconds = elapsed_time % 60;
    debug!("Running time: {hours}h:{minutes}m:{seconds}s");
}

/// Uses the clipped reads and vcf + json data if available to identify the coordinates of genomic
/// breakend locations, plus any connections between them that can be found via VCF.
fn identify_breakends_from_inputs(
    vcf_filename_opt: Option<PathBuf>,
    json_filename_opt: Option<PathBuf>,
    sample_id: String,
    exclude_regions: ExcludeRegions,
    target_opt: &Option<TargetCoordinate>,
    clipped_reads: HashMap<String, Vec<FwdStrandSplitReadSegment>>,
    allow_unphased: bool,
) -> (
    HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>,
    Vec<Connection>,
) {
    let vcf_coord_map: HashMap<String, Vec<Coordinate>>;
    let clip_coordinates: HashMap<Coordinate, HashSet<FwdStrandSplitReadSegment>>;
    let mut vcf_connections: Vec<Connection> = Vec::new();

    if let Some(vcf_filename) = vcf_filename_opt {
        (vcf_connections, vcf_coord_map) =
            get_vcf_breaks(vcf_filename, &exclude_regions, target_opt);
        let variant_ids: std::collections::HashSet<String> =
            vcf_coord_map.keys().cloned().collect();
        let mut variant_readnames_opt: Option<HashMap<String, Vec<String>>> = None;
        if let Some(json_filename) = json_filename_opt {
            variant_readnames_opt = Some(get_read_info_from_json(
                json_filename,
                sample_id,
                variant_ids,
            ));
        }
        clip_coordinates = cluster_support_builder::assign_clipped_reads_to_clusters(
            &clipped_reads,
            &vcf_coord_map,
            &variant_readnames_opt,
        );
    } else {
        // Find genomic break locations using clustered groups of clipped reads
        clip_coordinates = cluster_finder::find_breaks(&clipped_reads, allow_unphased);
    };

    (clip_coordinates, vcf_connections)
}

fn main() {
    ///////////////////////////////////////////////////////////////////////////
    // Set up
    let (args, target_opt) = set_up();
    let start_time = SystemTime::now();
    let mut exclude_regions: ExcludeRegions = load_exclude_regions(args.exclude_regions_path);

    ///////////////////////////////////////////////////////////////////////////
    // Get read data
    let clipped_reads: HashMap<String, Vec<FwdStrandSplitReadSegment>>;
    if let Some(target_coord) = target_opt.clone() {
        clipped_reads = get_split_alignments_from_region(args.bam_filename.clone(), target_coord);
    } else {
        clipped_reads = get_split_alignments(args.bam_filename.clone(), &mut exclude_regions);
    }

    let sample_id: String = get_sample_from_bam(args.bam_filename.clone());

    let (clip_coordinates, vcf_connections) = identify_breakends_from_inputs(
        args.vcf_filename,
        args.json_filename,
        sample_id,
        exclude_regions,
        &target_opt,
        clipped_reads,
        args.allow_unphased,
    );

    // Connect genomic break locations using connections from the VCF, alignments, & phasing
    let break_connections: HashMap<Connection, Vec<FwdStrandSplitReadSegment>> =
        connect_clusters(&clip_coordinates, &vcf_connections, args.allow_unphased);

    ///////////////////////////////////////////////////////////////////////////
    // Join 1-to-1 connections into full connected event graphs
    let event_graphs: Vec<EventGraph> = build_event_graphs(&break_connections, &clip_coordinates);

    ///////////////////////////////////////////////////////////////////////////
    //Annotate event graphs with ordering info and directionality
    let annotated_graphs: ComplexSVCalls = annotate_graphs(&event_graphs, &clip_coordinates);

    // Filter/write results
    ///////////////////////////////////////////////////////////////////////////
    let mut filtered_annotated_graphs: ComplexSVCalls;
    if args.nofilter || target_opt.is_some() {
        filtered_annotated_graphs = annotated_graphs;
    } else {
        filtered_annotated_graphs =
            block_filter::apply_filters(&annotated_graphs, args.bam_filename, args.max_coverage);
    }
    filtered_annotated_graphs.event_graphs.sort_unstable();
    result_writer::write_results(
        filtered_annotated_graphs,
        args.outdir.trim_end_matches('/').to_string(),
        args.prefix,
        args.write_unzipped,
    );
    log_time(start_time);
}
