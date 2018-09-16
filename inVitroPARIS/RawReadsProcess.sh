

### 1. split fastq

bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/splitFastq.pl -l TTG:59_48_2::CAA:59_48_3::ACC:766_72_1::CCG:766_72_3::others:unmatched -b 5:3 -d split_fastq -s hiseq_barcode.stat -U ZIKA_Raw_vitro.fastq"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/splitFastq.pl -l GTT:59_48_1::GCG:766_72_2::others:unmatched2 -b 6:3 -d split_fastq2 -s hiseq_barcode.stat -U unmatched.fastq"


### 2. read collapse

bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 59_48_1.fastq -o 59_48_1.uniq.fastq"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 59_48_2.fastq -o 59_48_2.uniq.fastq"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 59_48_3.fastq -o 59_48_3.uniq.fastq"

bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 766_72_1.fastq -o 766_72_1.uniq.fastq"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 766_72_2.fastq -o 766_72_2.uniq.fastq"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/readCollapse.pl -U 766_72_3.fastq -o 766_72_3.uniq.fastq"

### 3. read trimming

bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 59_48_1.uniq.fastq -o 59_48_1.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 59_48_2.uniq.fastq -o 59_48_2.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 59_48_3.uniq.fastq -o 59_48_3.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"

bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 766_72_1.uniq.fastq -o 766_72_1.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 766_72_2.uniq.fastq -o 766_72_2.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"
bsub -q Z-BNODE "perl ~/lipan/icSHAPE/icSHAPE/scripts/trimming.pl -U 766_72_3.uniq.fastq -o 766_72_3.trimmed.fastq -a ../../adaptor/icSHAPEparisAdapter.fa -l 13 -m 30"

### 4. remove rRNA reads

function remove_rRNA()
{
    input=$1
    output=$2
    rRNA_base=$3

    bsub -q Z-ZQF -n 20 -e error -o log \
        "bowtie2 \
            -p 20\
            -k 1000 \
            --very-sensitive-local \
            --local \
            -x $rRNA_base \
            -q \
            $input \
            | samtools view \
            | awk '\$2==4{printf \"@%s\n%s\n+\n%s\n\",\$1,\$10,\$11}' \
            > $output"

    #echo $CMD
    #bsub -q Z-BNODE -n 20 -e error -o log $CMD
}

IN=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/180912/split_fastq
OUT=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/180912/split_fastq/mapping
INDEX_rRNA=/150T/zhangqf/lipan/SMART_SHAPE/ref_sequence/human_rRNA/human_rRNA

remove_rRNA $IN/59_48_1.trimmed.fastq $OUT/59_48_1_rRNA.fastq $INDEX_rRNA
remove_rRNA $IN/59_48_2.trimmed.fastq $OUT/59_48_2_rRNA.fastq $INDEX_rRNA
remove_rRNA $IN/59_48_3.trimmed.fastq $OUT/59_48_3_rRNA.fastq $INDEX_rRNA

remove_rRNA $IN/766_72_1.trimmed.fastq $OUT/766_72_1_rRNA.fastq $INDEX_rRNA
remove_rRNA $IN/766_72_2.trimmed.fastq $OUT/766_72_2_rRNA.fastq $INDEX_rRNA
remove_rRNA $IN/766_72_3.trimmed.fastq $OUT/766_72_3_rRNA.fastq $INDEX_rRNA

### 5. map to transcriptome and virus

function mapping_transcriptome()
{
    input=$1
    output=$2
    index=$3

    CMD="STAR \
        --runMode alignReads \
        --genomeDir $index \
        --readFilesIn $input \
        --outFileNamePrefix $output \
        --outSAMattributes All \
        --alignIntronMin 1 \
        --alignIntronMax 999999999 \
        --scoreGapNoncan -4 \
        --scoreGapATAC -4 \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --runThreadN 20 \
        --outSAMtype SAM \
        --outFilterMismatchNoverLmax 0.04 \
        --limitOutSJcollapsed 10000000 \
        --limitIObufferSize 280000000 \
        --chimOutType SeparateSAMold \
        --outReadsUnmapped Fastx"

    echo $CMD
    bsub -q Z-ZQF -n 20 -e error -o log $CMD

}

IN=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/180912/split_fastq/mapping
OUT=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/180912/split_fastq/mapping
INDEX_59=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/index/59/
INDEX_766=/Share/home/zhangqf8/lipan/virus/PARIS/raw_data/index/766/

mapping_transcriptome $IN/59_48_1_rRNA.fastq $OUT/59_48_1/ $INDEX_59
mapping_transcriptome $IN/59_48_2_rRNA.fastq $OUT/59_48_2/ $INDEX_59
mapping_transcriptome $IN/59_48_3_rRNA.fastq $OUT/59_48_3/ $INDEX_59

mapping_transcriptome $IN/766_72_1_rRNA.fastq $OUT/766_72_1/ $INDEX_766
mapping_transcriptome $IN/766_72_2_rRNA.fastq $OUT/766_72_2/ $INDEX_766
mapping_transcriptome $IN/766_72_3_rRNA.fastq $OUT/766_72_3/ $INDEX_766


### 6. count virus reads

sam_group_trim -in 59_48_1/Aligned.out.sam,59_48_1/Chimeric.out.sam -out 59_48_1/59_48_1.sam
sam_group_trim -in 59_48_2/Aligned.out.sam,59_48_2/Chimeric.out.sam -out 59_48_2/59_48_2.sam
sam_group_trim -in 59_48_3/Aligned.out.sam,59_48_3/Chimeric.out.sam -out 59_48_3/59_48_3.sam

sam_group_trim -in 766_72_1/Aligned.out.sam,766_72_1/Chimeric.out.sam -out 766_72_1/766_72_1.sam
sam_group_trim -in 766_72_2/Aligned.out.sam,766_72_2/Chimeric.out.sam -out 766_72_2/766_72_2.sam
sam_group_trim -in 766_72_3/Aligned.out.sam,766_72_3/Chimeric.out.sam -out 766_72_3/766_72_3.sam


sam2dg -in 59_48_1/59_48_1.sam -out 59_48_1/59_48_1.dg -out_sam 59_48_1/59_48_1.dg.sam -min_overhang 1 -min_armlen 15
sam2dg -in 59_48_2/59_48_2.sam -out 59_48_2/59_48_2.dg -out_sam 59_48_2/59_48_2.dg.sam -min_overhang 1 -min_armlen 15
sam2dg -in 59_48_3/59_48_3.sam -out 59_48_3/59_48_3.dg -out_sam 59_48_3/59_48_3.dg.sam -min_overhang 1 -min_armlen 15

sam2dg -in 766_72_1/766_72_1.sam -out 766_72_1/766_72_1.dg -out_sam 766_72_1/766_72_1.dg.sam -min_overhang 1 -min_armlen 15
sam2dg -in 766_72_2/766_72_2.sam -out 766_72_2/766_72_2.dg -out_sam 766_72_2/766_72_2.dg.sam -min_overhang 1 -min_armlen 15
sam2dg -in 766_72_3/766_72_3.sam -out 766_72_3/766_72_3.dg -out_sam 766_72_3/766_72_3.dg.sam -min_overhang 1 -min_armlen 15

grep "KU501215.1" 59_48_1/59_48_1.dg | grep -v "ENST" | grep -v "-" | wc -l
grep "KU501215.1" 59_48_2/59_48_2.dg | grep -v "ENST" | grep -v "-" | wc -l
grep "KU501215.1" 59_48_3/59_48_3.dg | grep -v "ENST" | grep -v "-" | wc -l

grep "AY632535.2" 766_72_1/766_72_1.dg | grep -v "ENST" | grep -v "-" | wc -l
grep "AY632535.2" 766_72_2/766_72_2.dg | grep -v "ENST" | grep -v "-" | wc -l
grep "AY632535.2" 766_72_3/766_72_3.dg | grep -v "ENST" | grep -v "-" | wc -l

### 7. calculate RPKM

perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 59_48_1/Aligned.out.sam -o 59_48_1/rpkm.txt
perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 59_48_2/Aligned.out.sam -o 59_48_2/rpkm.txt
perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 59_48_3/Aligned.out.sam -o 59_48_3/rpkm.txt

perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 766_72_1/Aligned.out.sam -o 766_72_1/rpkm.txt
perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 766_72_2/Aligned.out.sam -o 766_72_2/rpkm.txt
perl ~/lipan/icSHAPE/icSHAPE/scripts/estimateRPKM.pl -i 766_72_3/Aligned.out.sam -o 766_72_3/rpkm.txt

grep "KU501215.1" */rpkm.txt
grep "AY632535.2" */rpkm.txt

### 8. Compile sam files

grep "KU501215.1" 59_48_1/59_48_1.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view -h > 59_48_vitro.sam
grep "KU501215.1" 59_48_2/59_48_2.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view >> 59_48_vitro.sam
grep "KU501215.1" 59_48_3/59_48_3.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view >> 59_48_vitro.sam
samtools view -bh 59_48_vitro.sam | samtools sort > 59_48_vitro.bam

grep "AY632535.2" 766_72_1/766_72_1.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view -h > 766_72_vitro.sam
grep "AY632535.2" 766_72_2/766_72_2.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view >> 766_72_vitro.sam
grep "AY632535.2" 766_72_3/766_72_3.sam | awk '($6~/M.*N.*M/&&$2==0)||(substr($0,1,1)=="@"){print $0}' | samtools view -bh | samtools sort | samtools view >> 766_72_vitro.sam
samtools view -bh 766_72_vitro.sam | samtools sort > 766_72_vitro.bam






