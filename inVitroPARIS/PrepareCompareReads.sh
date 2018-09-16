
### 1. sample sam file

samtools view /Users/lee/Virus/Upload_Data/PRVABC59_48_vitro.bam | wc -l    # 22093
samtools view /Users/lee/Virus/Upload_Data/PRVABC59_72_vivo.bam | wc -l     # 113678

samtools view /Users/lee/Virus/Upload_Data/MR766_72_vitro.bam | wc -l    # 7673
samtools view /Users/lee/Virus/Upload_Data/MR766_72_vivo.bam | wc -l     # 98015


samtools view -h /Users/lee/Virus/Upload_Data/PRVABC59_48_vitro.bam > PRVABC59_48_vitro.sam
samtools view -h -s 0.19434719118914828 /Users/lee/Virus/Upload_Data/PRVABC59_72_vivo.bam > PRVABC59_72_vivo.sam

samtools view -h /Users/lee/Virus/Upload_Data/MR766_72_vitro.bam > MR766_72_vitro.sam
samtools view -h -s 0.07828393613222466 /Users/lee/Virus/Upload_Data/MR766_72_vivo.bam > MR766_72_vivo.sam

### 2. sam2dg

sam2dg -in PRVABC59_48_vitro.sam -out PRVABC59_48_vitro.dg -min_overhang 1
sam2dg -in PRVABC59_72_vivo.sam -out PRVABC59_72_vivo.dg -min_overhang 1

sam2dg -in MR766_72_vitro.sam -out MR766_72_vitro.dg -min_overhang 1
sam2dg -in MR766_72_vivo.sam -out MR766_72_vivo.dg -min_overhang 1

### 3. Call dg clusters

dg_cluster -in PRVABC59_48_vitro.dg -out /dev/null -uniqDG yes -out_tab PRVABC59_48_vitro.tmp.tab -min_support 5
dg_cluster -in PRVABC59_72_vivo.dg -out /dev/null -uniqDG yes -out_tab PRVABC59_72_vivo.tmp.tab -min_support 5

grep "^>" PRVABC59_48_vitro.tmp.tab > PRVABC59_48_vitro.tab
grep "^>" PRVABC59_72_vivo.tmp.tab > PRVABC59_72_vivo.tab

rm PRVABC59_48_vitro.tmp.tab
rm PRVABC59_72_vivo.tmp.tab


dg_cluster -in MR766_72_vitro.dg -out /dev/null -uniqDG yes -out_tab MR766_72_vitro.tmp.tab -min_support 5
dg_cluster -in MR766_72_vivo.dg -out /dev/null -uniqDG yes -out_tab MR766_72_vivo.tmp.tab -min_support 5

grep "^>" MR766_72_vitro.tmp.tab > MR766_72_vitro.tab
grep "^>" MR766_72_vivo.tmp.tab > MR766_72_vivo.tab

rm MR766_72_vitro.tmp.tab
rm MR766_72_vivo.tmp.tab


### 4. sam file to matrix

sam2matrix -in PRVABC59_72_vivo.sam -out PRVABC59_72_vivo.matrix -min_overhang 2 -chr KU501215.1
sam2matrix -in PRVABC59_48_vitro.sam -out PRVABC59_48_vitro.matrix -min_overhang 2 -chr KU501215.1

sam2matrix -in MR766_72_vivo.sam -out MR766_72_vivo.matrix -min_overhang 2 -chr AY632535.2
sam2matrix -in MR766_72_vitro.sam -out MR766_72_vitro.matrix -min_overhang 2 -chr AY632535.2



