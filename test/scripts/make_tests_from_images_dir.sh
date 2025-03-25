cargo build
svtopo=target/debug/svtopo

vcf=test/data/vcf/sawfish.vcf.gz
reads_json=test/data/json/supporting_reads.json.gz
hg38_exclude=cnv.excluded_regions.hg38.bed.gz
wget https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz -O cnv.excluded_regions.hg38.bed.gz -nc

set -e
bam=$1
input_dir=$2
output_dir=$3

if [ ! -f "${bam}" ];
then
    echo "First argument must be a valid bam"
    exit
fi
if [ ! -d "${input_dir}" ];
then
    echo "Second argument must be a valid directory (containing svtopo images)"
    exit
fi

if [ ! -d "${output_dir}" ];
then
    echo "Third argument must be a valid directory"
    exit
fi


for image in "$input_dir"/*png;
do
    echo $image 
    image_start_region=$(basename $image | cut -f 1 -d '.' | cut -f 2 -d '_' | tr '_' '\n')
    image_start_region_underscores=$(basename $image | cut -f 1 -d '.' | cut -f 2 -d '_' | tr '_' '\n' | tr '-' '_')
    new_bam_name=$output_dir"/$image_start_region_underscores"_vcf.bam
    basename $image | cut -f 1 -d '.' | cut -f 2- -d '_' | tr '_' '\n' | tr '-' '\t' | bedtools sort -i - > .tmp.bed
    samtools view -hb $bam -ML .tmp.bed > $new_bam_name
    samtools index $new_bam_name
    image_prefix=$output_dir"/test"
    json_out_prefix=$output_dir"/"$image_start_region_underscores
    
    echo $image_start_region
    echo $image_start_region_underscores


    
    cmd="$svtopo \
        --bam $new_bam_name \
        --variant-readnames $reads_json \
        --vcf $vcf \
        --exclude-regions $hg38_exclude \
        --out-prefix $json_out_prefix \
        --write-unzipped"
    $cmd

    cmd="svtopovz \
        --json $json_out_prefix.json \
        --out-prefix $image_prefix"
    $cmd
done