use bgzip::{BGZFWriter, Compression};
use std::collections::{HashMap, HashSet};
use std::path::PathBuf;

use crate::containers::{ComplexSVCalls, Coordinate};
use flate2::write::GzEncoder;
use log::{debug, error, info};
use std::fs::File;
use std::io::Write;
use std::io::{self, BufWriter};

/// Converts the filtered block groups to a pretty JSON string and writes them out.
/// If a json outfile path has been input, it will be created/overwritten. Otherwise, writes
/// in un-compressed form to standard out.
pub fn write_results(
    filtered_block_groups: ComplexSVCalls,
    outdir: String,
    prefix: String,
    write_unzipped: bool,
) {
    debug!(
        "{} block groups after filtering",
        filtered_block_groups.event_graphs.len()
    );
    let (bed_path, json_path) =
        generate_output_paths(outdir.clone(), prefix.clone(), write_unzipped);
    let json_string = block_groups_to_json(&filtered_block_groups);

    if let Err(error) = write_json(json_string, json_path.clone()) {
        error!("Error writing JSON result to outdir {}\n{}", outdir, error);
        std::process::exit(exitcode::IOERR);
    }

    if write_unzipped {
        if let Err(error) = write_unzipped_bed(&filtered_block_groups, bed_path.clone()) {
            error!("Error writing BED result to outdir {}\n{}", outdir, error);
            std::process::exit(exitcode::IOERR);
        }
    } else if let Err(error) = write_gzipped_bed(&filtered_block_groups, bed_path.clone()) {
        error!("Error writing BED result to outdir {}\n{}", outdir, error);
        std::process::exit(exitcode::IOERR);
    }
    info!("JSON written to {}", json_path,);
    info!("BED written to {}", bed_path,);
}

/// Convert block groups into a pretty json String
fn block_groups_to_json(filtered_block_groups: &ComplexSVCalls) -> String {
    match serde_json::to_string_pretty(&filtered_block_groups) {
        Ok(json) => json,
        Err(_) => {
            error!("Failed to write blocks as JSON");
            std::process::exit(exitcode::IOERR);
        }
    }
}

/// Write the json file output
fn write_json(json_string: String, json_name: String) -> std::io::Result<()> {
    let json_outfile = PathBuf::from(json_name);
    let file_handle = File::create(&json_outfile)?;

    if json_outfile.extension().and_then(|ext| ext.to_str()) == Some("gz") {
        let mut gzip_filehandle = GzEncoder::new(file_handle, flate2::Compression::default());
        gzip_filehandle.write_all(json_string.as_bytes())?;
    } else {
        let mut writer = io::BufWriter::new(file_handle);
        writer.write_all(json_string.as_bytes())?;
        writer.flush()?
    }

    Ok(())
}

/// Given an already-validated directory path with a filename prefix,
/// generate the json.gz and bed.gz filenames for output.
fn generate_output_paths(outdir: String, prefix: String, write_unzipped: bool) -> (String, String) {
    if write_unzipped {
        let bed_filepath = format!("{}/{}.bed", outdir.clone(), prefix.clone());
        let json_filepath = format!("{}/{}.json", outdir.clone(), prefix.clone());
        (bed_filepath, json_filepath)
    } else {
        let bed_filepath = format!("{}/{}.bed.gz", outdir.clone(), prefix.clone());
        let json_filepath = format!("{}/{}.json.gz", outdir.clone(), prefix.clone());
        (bed_filepath, json_filepath)
    }
}

/// Given the blocks and vcf variant IDs mapped to break locations,
/// write a BED.GZ file that associates each block region spanned
/// to the image name that will correspond to it (one column) and
/// any vcf variant IDs that correspond to that image (one column, csv)
fn write_unzipped_bed(
    filtered_block_groups: &ComplexSVCalls,
    bed_name: String,
) -> std::io::Result<()> {
    let bed_path = PathBuf::from(bed_name);
    let mut bedfile = File::create(bed_path.clone())?;

    let version = env!("CARGO_PKG_VERSION");
    let header = format!(
        "#SVTopo v{}\n#chrom\tpos\tend\timage_name\tvariant_ids",
        version
    );
    writeln!(bedfile, "{}", header)?;

    for complex_event in filtered_block_groups.event_graphs.iter() {
        let mut coords: Vec<Coordinate> = Vec::new();
        let mut variant_ids: HashSet<String> = HashSet::new();

        for block in complex_event.iter() {
            if !block.coverages.is_empty() {
                coords.push(block.region.clone());
                variant_ids.extend(block.region.variant_ids.iter().cloned());
            }
        }
        let mut variant_ids_vec: Vec<String> = variant_ids.into_iter().collect();
        variant_ids_vec.sort_unstable();
        let image_name = get_image_name(&coords);
        let var_ids_string = variant_ids_vec.join(",");
        for coord in coords {
            let entry = format!(
                "{}\t{}\t{}\t{}\t{}",
                coord.start_chrom, coord.start, coord.end, image_name, var_ids_string
            );
            writeln!(bedfile, "{}", entry)?;
        }
    }

    Ok(())
}

/// Given the blocks and vcf variant IDs mapped to break locations,
/// write a BED.GZ file that associates each block region spanned
/// to the image name that will correspond to it (one column) and
/// any vcf variant IDs that correspond to that image (one column, csv)
fn write_gzipped_bed(
    filtered_block_groups: &ComplexSVCalls,
    bed_name: String,
) -> std::io::Result<()> {
    let bed_path = PathBuf::from(bed_name);
    let bedfile = File::create(bed_path.clone())?;

    // Wrap the file in a buffered writer
    let mut buf_writer = BufWriter::new(bedfile);
    let mut writer = BGZFWriter::new(&mut buf_writer, Compression::default());

    let version = env!("CARGO_PKG_VERSION");
    let header = format!(
        "#SVTopo v{}\n#chrom\tpos\tend\timage_name\tvariant_ids\n",
        version
    );

    writer.write_all(header.as_bytes())?;

    for complex_event in filtered_block_groups.event_graphs.iter() {
        let mut coords: Vec<Coordinate> = Vec::new();
        let mut variant_ids: Vec<String> = Vec::new();

        for block in complex_event.iter() {
            if !block.coverages.is_empty() {
                coords.push(block.region.clone());
                variant_ids.append(&mut block.region.variant_ids.clone());
            }
        }
        let mut variant_ids_vec: Vec<String> = variant_ids.into_iter().collect();
        variant_ids_vec.sort_unstable();
        let image_name = get_image_name(&coords);
        let var_ids_string = variant_ids_vec.join(",");
        for coord in coords {
            let entry = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                coord.start_chrom, coord.start, coord.end, image_name, var_ids_string
            );
            writer.write_all(entry.as_bytes())?;
        }
    }

    Ok(())
}

/// Reproduces the logic in the SVTopoVz get_image_name() function
/// by generating the name for the image using the regions that are included.
fn get_image_name(coords: &Vec<Coordinate>) -> String {
    let mut image_regions = HashMap::new();
    for coord in coords {
        image_regions
            .entry(coord.start_chrom.clone())
            .or_insert_with(Vec::new)
            .push(coord.start);
        image_regions
            .entry(coord.start_chrom.clone())
            .or_insert_with(Vec::new)
            .push(coord.end);
    }
    let mut chroms: Vec<&String> = image_regions.keys().collect();
    chroms.sort();
    let mut simplified_regions: Vec<String> = Vec::new();
    for chrom in chroms {
        if let Some(positions) = image_regions.get(chrom) {
            let start_opt = positions.iter().min();
            let end_opt = positions.iter().max();
            if let (Some(start), Some(end)) = (start_opt, end_opt) {
                let region_string = format!("{}-{}-{}", chrom, start, end);
                simplified_regions.push(region_string);
            }
        }
    }
    simplified_regions.join("__")
}
