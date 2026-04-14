########## Variant Calling  ########

### SNVs, clair3-nova  ###

SAMPLE0="$probandID"
SAMPLE1="$fatherID"
SAMPLE2="$motherID"

INPUT_DIR="/home/user" 
REF=${INPUT_DIR}/path/to/hg38/Homo_sapiens_assembly38.fasta              # change your reference file name here
OUTPUT_DIR="/path/to/snv-output/clair3-nova/${SAMPLE0}"        
THREADS="10"               
MODEL_C3="r1041_e82_400bps_hac_v430"         	  
MODEL_C3D="r1041_e82_400bps_hac_nova"      

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3-nova:latest \
  /opt/bin/run_clair3_nova.sh \
  --ref_fn=${INPUT_DIR}/path/to/hg38/Homo_sapiens_assembly38.fasta \
  --bam_fn_c=${INPUT_DIR}/path/to/SAMPLE0.bam \
  --bam_fn_p1=${INPUT_DIR}/path/to/SAMPLE1.bam \
  --bam_fn_p2=${INPUT_DIR}/path/to/SAMPLE2.bam \
  --sample_name_c=${SAMPLE0} \
  --sample_name_p1=${SAMPLE1} \
  --sample_name_p2=${SAMPLE2} \
  --threads=${THREADS} \
  --gvcf \
  --model_path_clair3="/opt/models/clair3_models/${MODEL_C3}" \
  --model_path_clair3_nova="/opt/models/clair3_nova_models/${MODEL_C3D}" \
  --output=${OUTPUT_DIR}      




bcftools view -i 'ALT!="." && FILTER="PASS"' SAMPLE0.vcf.gz -O z -o SAMPLE0_filtered.vcf.gz
bcftools view -i 'ALT!="." && FILTER="PASS"' SAMPLE1.vcf.gz -O z -o SAMPLE1_filtered.vcf.gz
bcftools view -i 'ALT!="." && FILTER="PASS"' SAMPLE2.vcf.gz -O z -o SAMPLE2_filtered.vcf.gz

############### UPLOAD SAMPLE0_filtered.vcf.gz SAMPLE1_filtered.vcf.gz SAMPLE2_filtered.vcf.gz on emedgene for tertiary analysis ################

# input: clair3-nova's output of ${SAMPLE[0]}.vcf.gz, ${SAMPLE[1]}.vcf.gz, ${SAMPLE[2]}.vcf.gz files
# output: merged vcf and de novo variants

# requires trio's ped file, reference sdf file
# example input
tee "${SAMPLE0}.ped" > /dev/null << FAM_CONTENT    #PED format pedigree
#fam-id/ind-id/pat-id/mat-id: 0=unknown
#sex: 1=male; 2=female; 0=unknown
#####change only the second last column of first row
#phenotype: -9=missing, 0=missing; 1=unaffected; 2=affected
#fam-id ind-id pat-id mat-id sex phen
${SAMPLE0} ${SAMPLE0} ${SAMPLE1} ${SAMPLE2} 2 0
${SAMPLE0} ${SAMPLE1} 0 0 1 1
${SAMPLE0} ${SAMPLE2} 0 0 2 1
FAM_CONTENT

# your reference sdf file path
REF_SDF_FILE_PATH=~/NAS_Nigro/ONT/MNESYS/MNESYS12_280725/clair3-nova/GRCh38.sdf

# output files
# merged and de novo vcfs
Merged_VCF=${SAMPLE0}_trio.vcf.gz
Merged_VCF_annotated=${SAMPLE0}_trio_ann.vcf.gz
denovo_VCF=${SAMPLE0}_trio_all_denovo.vcf.gz
denovo_VCF_sf=${SAMPLE0}_trio_high_quality_denovo.vcf.gz

# merge trio vcfs
bcftools merge ${OUTPUT_DIR}/${SAMPLE0}.vcf.gz \
${OUTPUT_DIR}/${SAMPLE1}.vcf.gz \
${OUTPUT_DIR}/${SAMPLE2}.vcf.gz \
--threads 32 -f PASS -0 -m all| bcftools view -O z -o ${OUTPUT_DIR}/${Merged_VCF}

# index
bcftools index ${OUTPUT_DIR}/${Merged_VCF}
#${BCFTOOLS} view ${M_VCF} -H | wc -l

# annotate with Mendelian inherrtance pattern
rtg mendelian -i ${OUTPUT_DIR}/${Merged_VCF} -o ${OUTPUT_DIR}/${Merged_VCF_annotated} --pedigree "${SAMPLE0}.ped" -t ${REF_SDF_FILE_PATH} |& tee MDL.log

# get de novo variants
bcftools view -i 'INFO/MCV ~ "0/0+0/0->0/1"' ${OUTPUT_DIR}/${Merged_VCF_annotated} -O z -o ${OUTPUT_DIR}/${denovo_VCF}
bcftools index ${OUTPUT_DIR}/${denovo_VCF}
# get high quality de novo variants
bcftools view -i "INFO/DNP>0.85" ${OUTPUT_DIR}/${denovo_VCF} -O z -o ${OUTPUT_DIR}/${denovo_VCF_sf}
bcftools index ${OUTPUT_DIR}/${denovo_VCF_sf}
# high quality de novo variants set is in ${denovo_VCF_sf}


