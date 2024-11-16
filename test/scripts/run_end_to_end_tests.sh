set -e
svtopo=$1
RED='\033[0;31m'
GREEN='\033[0;32m'
final_color=$GREEN


if [ "$svtopo" = "" ];
then
    echo "Missing svtopo binary"
    exit 1
fi

if [ ! -e "$svtopo" ];
then
    echo "Invalid path $svtopo"
    exit 1
fi

vcf=test/data/vcf/sawfish.vcf.gz
reads_json=test/data/json/supporting_reads.json.gz
wget https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz -O cnv.excluded_regions.hg38.bed.gz
hg38=cnv.excluded_regions.hg38.bed.gz

for bam in test/data/bam/*bam;
do
    echo $bam
    echo ""

    chrom=$(basename $bam | cut -f 1 -d '_')
    start=$(basename $bam | cut -f 2 -d '_')
    end=$(basename $bam | cut -f 3 -d '_' | cut -f 1 -d '.')
    runtype=$(basename $bam | cut -f 4 -d '_' | cut -f 1 -d '.')

    expected_image=$(ls test/data/image/"$chrom"-"$start"-"$end"*.png)
    expected_json="test/data/json/"$chrom"_"$start"_"$end".json"

    # run svtopo
    tmpname="test"
    test_json_name=${tmpname}_${chrom}_${start}.json
    if [ "$runtype" = "vcf" ];
    then
        $svtopo \
            --bam $bam \
            --verbose \
            --variant-readnames $reads_json \
            --vcf $vcf \
            --exclude-regions $hg38 \
            > $test_json_name \
            2> svtopo_tests.log
    else
        $svtopo \
            --bam $bam \
            --verbose \
            --exclude-regions $hg38 \
            > $test_json_name \
            2> svtopo_tests.log
    fi
    svtopovz \
        --json $test_json_name \
        --out-prefix $tmpname\
        --include-simple-dups\
        --include-simple-dels\
        --verbose 2>> svtopo_tests.log

    if [ $? -ne 0 ]; then 
        final_color="$RED"
    fi
    
    # compare results against saved files
    python test/scripts/compare_files.py -r $expected_json -t $test_json_name
    if [ $? -ne 0 ]; then 
        final_color="$RED"
    fi

    #test image output from svtopovz
    for test_image in ${tmpname}*.png;
    do
        test_image=$(ls $test_image)
        if [ $? -ne 0 ]; then 
            final_color="$RED"
        fi
        python test/scripts/compare_files.py -r $expected_image -t $test_image
        if [ $? -ne 0 ]; then 
            final_color="$RED"
        fi
        echo ""
    done
    rm ${tmpname}*png
    rm $test_json_name
    rm svtopo_tests.log
done

printf "\n${final_color}Tests finished\n"
