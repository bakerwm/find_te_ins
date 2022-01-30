#!/usr/bin/bash 

################################################################################
# extract TE insertions using long reads, with minimap2
#
################################################################################





################################################################################
# 1.alignment: minimap2
function align() {
    # input: fq.gz out_dir index
    # output: bam
    local fq=$1
    local out_dir=$2 # results/01.align
    local ref_fa=$3 # reference, fasta format
    local ref_name=$(basename ${ref_fa/.fa})
    local fname=$(basename ${fq/.fa})
    local fname=${fname/.fq.gz}
    local bam="${out_dir}/${fname}.bam"
    local log="${out_dir}/run_minimap2.${ref_name}.stderr"
    # local ref_fa="/data/yulab/wangming/data/genome/dm6/bigZips/dm6.fa"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    [[ ! -f ${bam} ]] && minimap2 -t 12 -ax map-ont ${ref_fa} ${fq} 2> ${log} | samtools view -bhS - | samtools sort -o ${bam} - && samtools index ${bam}
    echo $bam
}


################################################################################
# 2.extract raw insertions: from CIGAR
function get_raw_ins() {
    # input: bam out_dir
    # output: ins.bed
    local bam=$1
    local out_dir=$2 # results/02.raw_ins
    local fname=$(basename ${bam/.bam})
    local out_bed="${out_dir}/${fname}.raw_ins.bed"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # run
    python scripts/extract_INS/extract_ins.py ${bam} > ${out_bed}
    echo ${out_bed}
}


################################################################################
# 04.annotation of ins (te names)
# python scripts/extract_INS/anno_te.py -x dm6_transposon.fa.fai ${bam} | sort -k4,4 -k5,5n > ${out}
# python scripts/extract_INS/merge_te.py -w 100 $t > $out

################################################################################
# get insertion, full version

function get_ins() {
    local fq=$1
    local out_dir=$2
    local window=100 #
    local src_dir="/data/yulab/wangming/work/yu_2021/piwi_lxh/results/20220123_ONT_TE_insertion/scripts/extract_INS"
    local fname=$(basename ${fq/.fq.gz})
    local prj_dir="${out_dir}/${fname}"
    [[ ! -d ${prj_dir} ]] && mkdir -p ${prj_dir}
    
    ## 1. align to reference genome
    echo "[1/8] align to reference genome"
    local ref_fa="/data/yulab/wangming/data/genome/dm6/bigZips/dm6.fa"
    bam=$(align ${fq} ${prj_dir} ${ref_fa})
    
    ## 2. extract raw insertions
    echo "[1/8] extract raw insertions from BAM, by CIGAR"
    local raw_ins_bed="${prj_dir}/${fname}.raw_ins.bed"
    [[ ! -f ${raw_ins_bed} ]] && python ${src_dir}/extract_ins.py ${bam} > ${out_bed}
    
    ## 3. convert insertion to fasta
    echo "[1/8] convert raw insertions to fasta"
    local raw_ins_fa="${prj_dir}/${fname}.raw_ins.fa"
    [[ ! -f ${raw_ins_fa} ]] && python ${src_dir}/bed2fa.py ${raw_ins_bed}
    
    ## 4. align insertions to transposon
    echo "[1/8] align raw_insertion to transposon"
    local te_fa="/data/yulab/wangming/data/genome/dm6/dm6_transposon/dm6_transposon.fa"
    te_bam=$(align ${raw_ins_fa} ${prj_dir} ${te_fa})
    
    ## 5. add te name to insertion
    echo "[1/8] extract transposon name for insertions"
    local te_fa_fai=${te_fa}.fai
    local te_ins_txt="${prj_dir}/${fname}.te_ins.raw.txt"
    [[ ! -f ${te_ins_txt} ]] && python ${src_dir}/anno_te.py -x ${te_fa_fai} ${te_bam} | sort -k4,4 -k5,5n > ${te_ins_txt}
    
    ## 6. merge te insertions by window
    echo "[1/8] merge raw_insertions by window=${window}"
    local te_ins_bed="${prj_dir}/${fname}.te_ins.bed"
    python ${src_dir}/merge_te.py -w ${window} ${te_ins_txt} > ${te_ins_bed}

    ## 7. quant te insertions
    echo "[1/8] get the read count for each insertion: ~5min"
    local te_ins_gtf="${prj_dir}/${fname}.te_ins.gtf"
    local te_ins_quant="${prj_dir}/${fname}.te_ins.quant.txt"
    local quant_stdout="${prj_dir}/${fname}.te_ins.quant.stdout"
    local quant_stderr="${prj_dir}/${fname}.te_ins.quant.stderr"
    [[ ! -f ${te_ins_gtf} ]] && python ${src_dir}/bed2gtf.py ${te_ins_bed} > ${te_ins_gtf}
    [[ ! -f ${te_ins_quant} ]] && featureCounts -L -a ${te_ins_gtf} -o ${te_ins_quant} -t gene ${bam} 1>${quant_stdout} 2> ${quant_stderr}
    
    ## 8. add quant to te insertions
    echo "[1/8] write find insertions to file"
    local te_ins_final="${prj_dir}/${fname}.te_ins.final.bed"
    [[ ! -f te_ins_final ]] && python ${src_dir}/quant_bed.py ${te_ins_bed} ${te_ins_quant} > ${te_ins_final}

    ## 9. Done
    echo "Done!"
}


## run
for fq in data/raw_data/ONT*gz
do
    get_ins ${fq} results/extract_INS/
#     break
done


