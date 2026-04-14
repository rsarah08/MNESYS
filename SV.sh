SAMPLE0="$probandID"
SAMPLE1="$fatherID"
SAMPLE2="$motherID"

BAM_PROBAND="/path/to/SAMPLE0.bam"
BAM_FATHER="/path/to/SAMPLE1.bam"
BAM_MOTHER="/path/to/SAMPLE2.bam"

#sniffles --reference /path/to/hg38/Homo_sapiens_assembly38.fasta --input ${BAM_PROBAND} --tandem-repeats /path/to/Sniffles/annotations/human_GRCh38_no_alt_analysis_set.trf.bed --snf ${SAMPLE0}.SNIFFLES2.snf --vcf ${SAMPLE0}.SNIFFLES2.vcf --sample-id ${SAMPLE0}
#sniffles --reference /path/to/hg38/Homo_sapiens_assembly38.fasta --input ${BAM_FATHER} --tandem-repeats /path/to/Sniffles/annotations/human_GRCh38_no_alt_analysis_set.trf.bed --snf ${SAMPLE1}.SNIFFLES2.snf --vcf ${SAMPLE1}.SNIFFLES2.vcf --sample-id ${SAMPLE1}
#sniffles --reference /path/to/hg38/Homo_sapiens_assembly38.fasta --input ${BAM_MOTHER} --tandem-repeats /path/to/Sniffles/annotations/human_GRCh38_no_alt_analysis_set.trf.bed --snf ${SAMPLE2}.SNIFFLES2.snf --vcf ${SAMPLE2}.SNIFFLES2.vcf --sample-id ${SAMPLE2}

#sniffles --input  {SAMPLE0}.SNIFFLES2.snf ${SAMPLE1}.SNIFFLES2.snf ${SAMPLE2}.SNIFFLES2.snf --vcf ${SAMPLE0}_trio.vcf

#bgzip ${SAMPLE0}_trio.vcf
#tabix -p vcf ${SAMPLE0}_trio.vcf.gz

#bgzip ${SAMPLE0}.SNIFFLES2.vcf
#bgzip ${SAMPLE1}.SNIFFLES2.vcf
#bgzip ${SAMPLE2}.SNIFFLES2.vcf

#tabix -p vcf  ${SAMPLE0}.SNIFFLES2.vcf.gz
#tabix -p vcf  ${SAMPLE1}.SNIFFLES2.vcf.gz
#tabix -p vcf  ${SAMPLE2}.SNIFFLES2.vcf.gz

source $activate cutefc_env

cuteSV ${BAM_PROBAND} ~/path/to/hg38/Homo_sapiens_assembly38.fasta ${SAMPLE0}_cuteSV.vcf ${SAMPLE0}/ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads 32 --genotype --sample ${SAMPLE0} --max_size -1
#cuteSV ${BAM_FATHER} ~/path/to/hg38/Homo_sapiens_assembly38.fasta ${SAMPLE1}_cuteSV.vcf ${SAMPLE1}/ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads 32 --genotype --sample ${SAMPLE1} --max_size -1
cuteSV ${BAM_MOTHER} ~/path/to/hg38/Homo_sapiens_assembly38.fasta ${SAMPLE2}_cuteSV.vcf ${SAMPLE2}/ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --threads 32 --genotype --sample ${SAMPLE2} --max_size -1

source $activate nanopore 

bgzip ${SAMPLE0}_cuteSV.vcf
bgzip ${SAMPLE1}_cuteSV.vcf
bgzip ${SAMPLE2}_cuteSV.vcf

tabix -p vcf  ${SAMPLE0}_cuteSV.vcf.gz
tabix -p vcf  ${SAMPLE1}_cuteSV.vcf.gz
tabix -p vcf  ${SAMPLE2}_cuteSV.vcf.gz

bcftools merge ${SAMPLE0}_cuteSV.vcf.gz  ${SAMPLE1}_cuteSV.vcf.gz ${SAMPLE2}_cuteSV.vcf.gz -Oz -o  ${SAMPLE0}_trio_cuteSV.vcf.gz
tabix -p vcf ${SAMPLE0}_trio_cuteSV.vcf.gz
mkdir ${SAMPLE0}_isec_sniffles_cutesv
bcftools isec -c all ${SAMPLE0}_trio.vcf.gz ${SAMPLE0}_trio_cuteSV.vcf.gz -p ${SAMPLE0}_isec_sniffles_cutesv
cd ${SAMPLE0}_isec_sniffles_cutesv
sed -i 's/.SNIFFLES2//g' 0002.vcf
bgzip 0002.vcf
bgzip 0003.vcf
tabix -p vcf 0002.vcf.gz
tabix -p vcf 0003.vcf.gz
bcftools concat 0002.vcf.gz 0003.vcf.gz -Oz -o ${SAMPLE0}_trio_combined_concat.vcf.gz -a
bcftools norm -m -any -Oz -o ${SAMPLE0}_trio_combined_concat.split.vcf.gz ${SAMPLE0}_trio_combined_concat.vcf.gz
tabix -p vcf ${SAMPLE0}_trio_combined_concat.split.vcf.gz
truvari collapse --write-resolved --reference /path/to/hg38/Homo_sapiens_assembly38.fasta -i ${SAMPLE0}_trio_combined_concat.split.vcf.gz -o ${SAMPLE0}_trio_combined_concat.split.truvari.vcf.gz -c ${SAMPLE0}_trio_combined_truvari_collapsed.vcf --sizemax -1 --sizemin 30 --refdist 500 --pctseq 0.5 --pctsize 0.2 --pctovl 0.2
kanpig trio --input ${SAMPLE0}_trio_combined_concat.split.truvari.vcf.gz --proband ${BAM_PROBAND} --father ${BAM_FATHER} --mother ${BAM_MOTHER} --reference /path/to/hg38/Homo_sapiens_assembly38.fasta --proband-sample ${SAMPLE0} --father-sample ${SAMPLE1} --mother-sample ${SAMPLE2} --out ${SAMPLE0}_trio_combined_concat_split_truvari_kanpig.vcf
bgzip ${SAMPLE0}_trio_combined_concat_split_truvari_kanpig.vcf
tabix -p vcf ${SAMPLE0}_trio_combined_concat_split_truvari_kanpig.vcf.gz

for file in *truvari_kanpig.vcf.gz; do     echo "Processing $file..."
bcftools sort "$file" -O z -o "${file%.vcf}.truvari.kanpig.sorted.vcf.gz";
tabix -p vcf "${file%.vcf}.truvari.kanpig.sorted.vcf.gz"; done

##################### ANNOTATION #############################


ls -d "$PWD"/${SAMPLE0}_trio_combined_concat_split_truvari_kanpig.vcf.gz.truvari.kanpig.sorted.vcf.gz >needLR.txt
mkdir needLR_output
~/path/to/needLR_local/needLR_3.4_basic.sh -g /path/to/hg38/Homo_sapiens_assembly38.fasta -l needLR.txt

AnnotSV     -annotationMode both     -SVinputFile ${SAMPLE0}_trio_combined_concat_split_truvari_kanpig.vcf.gz.truvari.kanpig.sorted.vcf.gz     -genomeBuild GRCh38     -hpo HP:0007302     -outputDir ${SAMPLE0}_trio_isec_annotsv     -outputFile ${SAMPLE0}_trio_combined_concat_truvari_kanpig_annotsv     -SVminSize 0     -annotationsDir /path/toAnnotSV/share/AnnotSV     -variantconvertDir /path/tovariantconvert     -vcf 1     -variantconvertMode combined     -SVinputInfo 1

