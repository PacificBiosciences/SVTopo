set -e
cargo build --release
cd SVTopoVz/
python setup.py install
cd ..
pytest
svtopo=target/release/svtopo
RED='\033[0;31m'
GREEN='\033[0;32m'
final_color=$GREEN

vcf=test/data/vcf/sawfish.vcf.gz
reads_json=test/data/json/supporting_reads.json.gz
hg38_exclude=cnv.excluded_regions.hg38.bed.gz
wget https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.hg38.bed.gz -O cnv.excluded_regions.hg38.bed.gz

for bam in test/data/bam/*bam;
do
    echo $bam
    echo ""

    tmp_test_dir=".tmp_test"
    prefix="tmp"
    chrom=$(basename $bam | cut -f 1 -d '_')
    start=$(basename $bam | cut -f 2 -d '_')
    end=$(basename $bam | cut -f 3 -d '_' | cut -f 1 -d '.')
    runtype=$(basename $bam | cut -f 4 -d '_' | cut -f 1 -d '.')
    if [ "$start" -eq 189119276 ] || [ "$start" -eq 87906696 ] || [ "$start" -eq 65512263 ] || [ "$start" -eq 174937531 ]; 
    then
        continue
    fi

    # if [ ! "$start" -eq 138232624 ]; 
    # then
    #     continue
    # fi
    expected_image=$(ls test/data/image/"$chrom"-"$start"-"$end"*.png)
    expected_html="test/data/site_pages/"$chrom"-"$start"-"$end".html"
    expected_json="test/data/json/"$chrom"_"$start"_"$end".json"
    expected_bed="test/data/bed/"$chrom"_"$start"_"$end".bed"
    rm -rf $tmp_test_dir
    mkdir -p $tmp_test_dir

    # run svtopo
    tmpname="$tmp_test_dir/$prefix"
    test_html="$tmp_test_dir/index.html"
    test_json_name=${tmpname}.json
    test_bed_name=${tmpname}.bed
    if [ "$runtype" = "vcf" ];
    then
        cmd="$svtopo \
            --bam $bam \
            --variant-readnames $reads_json \
            --vcf $vcf \
            --exclude-regions $hg38_exclude \
            --prefix $prefix \
            --svtopo-dir $tmp_test_dir \
            --write-unzipped"
        echo "$cmd"
        $cmd
    else
        cmd="$svtopo \
            --bam $bam \
            --exclude-regions $hg38_exclude \
            --prefix $prefix \
            --svtopo-dir $tmp_test_dir 
            --write-unzipped"
        echo "$cmd"
        $cmd
    fi
    cmd="svtopovz \
        --svtopo-dir $tmp_test_dir" \

    echo "$cmd"
    $cmd

    if [ $? -ne 0 ]; then 
        final_color="$RED"
    fi
    
    # compare results against saved files
    python test/scripts/compare_files.py -r $expected_json -t $test_json_name
    if [ $? -ne 0 ]; then 
        final_color="$RED"
    fi

    python test/scripts/compare_files.py -r $expected_bed -t $test_bed_name
    if [ $? -ne 0 ]; then 
        final_color="$RED"
    fi
    
    python test/scripts/compare_files.py -r $expected_html -t $test_html
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
    # only to be used for renewing tests after a change has been validated
    # cp $test_json_name $expected_json
    # cp $test_bed_name $expected_bed
    # cp $test_image $expected_image
    # cp $test_html $expected_html
    rm -rf $tmp_test_dir
done

printf "\n${final_color}Tests finished\n"