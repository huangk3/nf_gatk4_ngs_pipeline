#!/bin/env nextflow

if (! params.keySet().containsAll(['batch'])) {
  println "please provide the parameter 'batch'"
  println 'nextflow run -c config_file pipeline.nf --batch batchNum'
  exit 1
}


REF = file(params.ref)
REF_DIC = file(params.ref_dic)
AXIOM = file(params.axiom)
MILLS = file(params.mills)
OneThousand = file(params.one_thousand)
OMINI = file(params.omini)
HAPMAP = file(params.hapmap)
TARGETS = file(params.targets)
BAITS = file(params.baits)
CONTROL_EXOMES = file(params.controlExomes).readLines()
OUTDIR = file(params.output_dir)
DBSNP = file(params.dbsnp)
snp_recal_tranche_values = params.snv_recal_tranche_values
snp_recal_anno_values = params.snv_recal_annotation_values
indel_recal_tranche_values = params.indel_recal_tranche_values
indel_recal_anno_values = params.indel_recal_annotation_values
BATCH = params.batch
dbNSFP_dir = file("/sysapps/cluster/software/VEP/84-goolf-1.7.20/.vep/dbNSFP/2.9.3")

logParams(params, "nextflow_parameters.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GRIS WES pipeline for Exome-Seq Preprocessing v${VERSION}"
log.info "Nextflow Version:     $workflow.nextflow.version"
log.info "Command Line:         $workflow.commandLine"
log.info "Batch:                ${BATCH}"
log.info "========================================="

//Channel.fromPath("/nethome/huangk3/bcbb/Andrew_Demo/WESTrainingCourse/daughter*/*{1,2}.fq").into{bams; extractSampleSizes}
Channel.fromPath("/nethome/huangk3/bcbb/DATA_ANALYSIS/nextflow/BAMs/P00*.bam").into{bams; extractSampleSizes}


N=extractSampleSizes.count()

process runBam2Fastq {
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/Bam2Fastq/${sampleID}"
  conda "samtools"
  input:
    file(bam) from bams
  output:
    tuple sampleID, file(fq1), file(fq2), file(unpaired_reads) into runBam2FastqOut
  script:
    sampleID = bam.getSimpleName()
    fq1 = sampleID + ".R1.fq.gz"
    fq2 = sampleID + ".R2.fq.gz"
    unpaired_reads = sampleID + ".singleton.fq.gz"
      """
        samtools fastq -t -@4 -1 ${fq1} -2 ${fq2} -s ${unpaired_reads} ${bam}
      """
}

// ------------------------------------------------------------------------------------------------------------
//
// Preprocess reads
// 1) Align fastq to reference
// 2) sort and index the bam
//
// ------------------------------------------------------------------------------------------------------------



process runFastqToBam {
  tag "${BATCH}|${sampleID}"
  label "highcpu"
  publishDir "${OUTDIR}/${BATCH}/FastqToSam/${sampleID}/"
  conda "bwa picard"
  input: 
    tuple sampleID, file(fq1), file(fq2), file(singleton_reads) from runBam2FastqOut
  output:
    tuple sampleID, file(outBam), file(outBai) into (runFastqToBamOut, runFastQCIn)
    val sampleID into sampleID_list, sampleID_mapping
  script:
    outBam = sampleID + ".sorted.bam"
    outBai = sampleID + ".sorted.bai"
    """
      bwa mem -t 8 -R "@RG\\tID:${sampleID}\\tLB:${sampleID}\\tSM:${sampleID}\\tPL:ILLUMINA" ${REF}  ${fq1} ${fq2} | \
      picard -Xmx32G SortSam I=/dev/stdin O=${outBam} \
      MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate
    """
}

process runFastQC {
  tag "${BATCH}|${sampleID}"
  publishDir "${OUTDIR}/${BATCH}/FastQC/${sampleID}/"
  conda "fastqc"	    
  input:
    tuple sampleID, file(bam), file(bai) from runFastQCIn

  output:
    tuple sampleID, file("*.zip") into FastQCOut
    file("*.html") into FastQCOut2
   	
  script:

    """
      fastqc -t 1 -o . -f bam ${bam}
    """
}


process runMarkDuplicates {
  tag "${BATCH}|${sampleID}"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/MarkDuplicate/${sampleID}/"
  conda "picard"
  input:
    tuple sampleID, file(bam), file(bai) from runFastqToBamOut

  output:
    tuple sampleID, file(outBam) , file(outBai), file(outMetrics) into (runMarkDuplicateOut, runDepthOfCoverageIn)
    tuple BATCH, file(outBam) into runFastQC
  script:
    outBam = sampleID + ".markdup.bam"
    outBai = sampleID + ".markdup.bai"
    outMetrics = sampleID + ".markdup.metrics.txt" 
  """
    picard -Xmx32G MarkDuplicates I=${bam} O=${outBam} \
    VALIDATION_STRINGENCY=LENIENT AS=true CREATE_INDEX=true M=${outMetrics} 
  """
}


process runDepthOfCoverage {
  tag "${BATCH}|${sampleID}"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/DepthOfCoverage/${sampleID}/"
  input:
    tuple sampleID, file(bam), file(bai), file(outMetrics) from runDepthOfCoverageIn

  output:
    file("${prefix}*") into runDepthOfCoverageOut
  script:
    prefix = sampleID + "."
  """
    module purge
    module load gatk/3.8.1-Java-1.8.0_92
    java -Xmx32G -jar \${EBROOTGATK}/GenomeAnalysisTK.jar -R ${REF} -T DepthOfCoverage -I ${bam} -geneList:REFSEQ ${params.geneList} -mmq 20 -mbq 10 -omitBaseOutput \
    --outputFormat csv -L ${TARGETS} -ct 8 -ct 10 -ct 20 -ct 100 -o ${sampleID}
  """
}


process runBaseRecalibrator {
  tag "$BATCH|$sampleID"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/BQSR/${sampleID}/"
  conda "gatk4"
  input:
    tuple sampleID, file(bam) , file(bai), file(mkdup_metrics) from runMarkDuplicateOut

  output:
    tuple sampleID, file(bam) , file(bai), file(recal_table) into runBaseRecalibratorOut

  script:
    recal_table=sampleID + ".recal_data.table"
    
    """
      gatk -Xmx32G BaseRecalibrator -R ${REF} -I ${bam} -O ${recal_table} --known-sites ${DBSNP} --known-sites ${OneThousand} -known-sites ${MILLS}
    """
}

process runApplyBQSR {
  tag "$BATCH|$sampleID"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/BQSR/${sampleID}/"
  conda "gatk4"
  input:
    tuple sampleID, file(bam) , file(bai), file(recal_table) from runBaseRecalibratorOut

  output:
    tuple sampleID, file(outBam), file(outBai) into runApplyBQSROut
   
  script:
    outBam = sampleID + ".bqsr.bam"
    outBai = sampleID + ".bqsr.bai"

    """
      gatk -Xmx32G ApplyBQSR -R ${REF} -I ${bam} -bqsr ${recal_table} -O ${outBam}
    """
}

Channel.from('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY').into{chr_gt; chr_vep}

process runHaplotypeCaller {
  tag "$BATCH|$sampleID"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/gVCF/${sampleID}/"
  conda "gatk4"
  input:
    tuple sampleID, file(bam) , file(bai) from runApplyBQSROut

  output:
    file(gvcf) into runHaplotypeCallerVcfOut
    file(gvcf_index) into runHaplotypeCallerVcfIndexOut

  script:
    gvcf = sampleID + ".g.vcf.gz"
    gvcf_index = sampleID + ".g.vcf.gz.tbi"

    """
      gatk -Xmx32G HaplotypeCaller -R ${REF} -I ${bam} -O ${gvcf} -stand-call-conf 30 --dbsnp ${DBSNP} -L ${TARGETS} -ERC GVCF
      tabix -p vcf -f ${gvcf}
    """
}

sampleID_list.collectFile(newLine: true, storeDir: "${OUTDIR}/${BATCH}") {item->
  ["${BATCH}.sampleList.txt", item]
}

//sampleID_mapping.collectFile(newLine: true, storeDir: "${OUTDIR}/${BATCH}") {item->
//  ["${BATCH}.gVCF.mapping.txt", item + "\t" + "${OUTDIR}/${BATCH}/"+ item +"/HaplotypeCaller/" + item + ".g.vcf.gz"]
//}
//--sample-name-map ${OUTDIR}/${BATCH}/${BATCH}.gVCF.mapping.txt

//runHaplotypeCallerOut_grouped_by_batch = runHaplotypeCallerOut.groupTuple(by: [0])
//total_intervals = Channel.fromPath("${TARGETS}").countLines()
//process runDynamicallyCombineIntervals {
//  publishDir "${OUTDIR}/${BATCH}/splittedIntervals/"
//  input:
//    val n from N
//    val total from total_intervals

//  output:
//    file("p*.intervals") into runDynamicallyCombineIntervalsOut

//  script:
//    lines_per_file =  (1 + total / (n + 30)).toInteger()
//    """
//      split -d -a4 -l $lines_per_file --additional-suffix ".intervals" ${TARGETS} "p"
//    """
//}


process runGenomicsDBImport {
  tag "${BATCH}|${chr}"
  label "highcpu"
  publishDir "${OUTDIR}/${BATCH}/GenomicsDB"
  conda "gatk4"
  input:
    file(gvcf) from runHaplotypeCallerVcfOut.collect() 
    file(gvcf_index) from runHaplotypeCallerVcfIndexOut.collect()
    each chr from chr_gt
  output:
    tuple BATCH, chr, file(genomeDb) into runGenomicsDBImportOut

  script:
//  p = interval.getSimpleName()
    genomeDb = BATCH + ".chr" + chr + ".genomeDb"
  """
    gatk -Xmx32G GenomicsDBImport --genomicsdb-workspace-path ${genomeDb} --overwrite-existing-genomicsdb-workspace true \
    --batch-size 50 -L ${chr} ${gvcf.collect { "-V $it " }.join()}
  """ 
}

process runGenotypeGVCFs {
  tag "$BATCH|${chr}"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/rawVcfByChr"
  conda "gatk4"
  input:
    tuple BATCH, chr, file(genomeDb) from runGenomicsDBImportOut

  output:
    tuple BATCH, chr, file(outVcf), file(outVcf_index) into runGenotypeGVCFsOut

  script:
    outVcf = BATCH + ".chr" + chr + ".vcf.gz"
    outVcf_index = BATCH + ".chr" + chr + ".vcf.gz.tbi"
  """
    gatk -Xmx32G GenotypeGVCFs -R ${REF} -V gendb://${genomeDb} -O ${outVcf} -G StandardAnnotation --only-output-calls-starting-in-intervals -L ${chr}
  """
}

runGenotypeGVCFsOut_group_by_batch=runGenotypeGVCFsOut.groupTuple(by: [0])

process runSortVcfs {
  tag "$BATCH"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  conda "gatk4"
  input:
    tuple BATCH, chr_list, vcf_list, vcf_index_list from runGenotypeGVCFsOut_group_by_batch
  output:
    tuple file(outVcf), file(outVcf_index) into runSortVcfsOut, ApplyVQSR  
  script:
    outVcf = BATCH + ".vcf.gz"
    outVcf_index = BATCH + ".vcf.gz.tbi"
    
  """
    gatk -Xmx32G SortVcf --INPUT ${vcf_list.join(" --INPUT ")} --OUTPUT ${outVcf} --SEQUENCE_DICTIONARY ${REF_DIC}
  """
} 

process runHardFilterAndMakeSitesOnlyVcf {
  tag "$BATCH"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  conda "gatk4 picard"
  input:
    tuple file(vcf), file(vcf_index) from runSortVcfsOut
  output:
    tuple file(sitesOnly_vcf), file(sitesOnly_vcf_index) into (snvRecal, indelRecal)
  script:
    sitesOnly_vcf = BATCH+".siteOnly.vcf.gz"
    sitesOnly_vcf_index = BATCH+".siteOnly.vcf.gz.tbi"
    filtered_vcf = BATCH+".filtered.vcf.gz"
    filtered_vcf_index = BATCH+".filtered.vcf.gz.tbi"
    """
      gatk -Xmx32G VariantFiltration --filter-expression "ExcessHet > ${params.excess_het_threshold}" --filter-name ExcessHet -O ${filtered_vcf} -V ${vcf}
      picard -Xmx32G MakeSitesOnlyVcf INPUT= ${filtered_vcf} OUTPUT=${sitesOnly_vcf}
    """
}

process runIndelRecal {
  tag "$BATCH"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/VQSR"
  conda "gatk4"
  input:
    tuple file(vcf), file(vcf_index) from indelRecal
  output:
    tuple file(indel_recal), file(indel_recal_index), file(indel_tranches) into runIndelRecalOut
  script:
    indel_recal = BATCH + ".indel.recal"
    indel_recal_index = BATCH + ".indel.recal.idx"
    indel_tranches = BATCH + ".indel.tranches"
  """
    module load gatk/4.0-8-Java-1.8.0_92
    \$EBROOTGATK/gatk -Xmx32G VariantRecalibrator -R ${REF} -V ${vcf} --output ${indel_recal} --tranches-file ${indel_tranches} --trust-all-polymorphic \
    -an ${indel_recal_anno_values.join(" -an ")} -tranche ${indel_recal_tranche_values.join(" -tranche ")} \
    -mode INDEL --max-gaussians 4 --resource mills,known=false,training=true,truth=true,prior=12:${MILLS} \
    --resource axiomPoly,known=false,training=true,truth=false,prior=10:${AXIOM} \
    --resource dbsnp,known=true,training=false,truth=false,prior=2:${DBSNP}
  """
}

process runSnvRecal {
  tag "$BATCH"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/VQSR"
  
  echo true
  input:
    tuple file(vcf), file(vcf_index) from snvRecal
  output:
    tuple file(snv_recal), file(snv_recal_index), file(snv_tranches) into runSnvRecalOut
  script:
    snv_recal = BATCH + ".snv.recal"
    snv_recal_index = BATCH + ".snv.recal.idx"
    snv_tranches = BATCH + ".snv.tranches"
  """
    module load gatk/4.0-8-Java-1.8.0_92
    \$EBROOTGATK/gatk -Xmx32G --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' VariantRecalibrator -R ${REF} -V ${vcf} --output ${snv_recal} --tranches-file ${snv_tranches} --trust-all-polymorphic \
    -an ${snp_recal_anno_values.join(" -an ")} -tranche ${snp_recal_tranche_values.join(" -tranche ")} \
    -mode SNP --max-gaussians 6 --resource hapmap,known=false,training=true,truth=true,prior=15:${HAPMAP} \
    --resource omni,known=false,training=true,truth=true,prior=12:${OMINI} \
    --resource 1000G,known=false,training=true,truth=false,prior=10:${OneThousand} \
    --resource dbsnp,known=true,training=false,truth=false,prior=7:${DBSNP}
  """
}


process runApplyVQSR {
  tag "$BATCH"
  label "highcpu"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  conda "gatk4"
  input:
    tuple file(indel_recal), file(indel_recal_index), file(indel_tranches) from runIndelRecalOut
    tuple file(snv_recal), file(snv_recal_index), file(snv_tranches) from runSnvRecalOut
    tuple file(vcf), file(vcf_index) from ApplyVQSR
  output:
    tuple file(recal_vcf), file(recal_vcf_index) into (runApplyRecalibrationOut, runNormalizationIn)
  script:
    recal_vcf = BATCH + ".recal.vcf.gz"
    recal_vcf_index = BATCH + ".recal.vcf.gz.tbi"
    """
      module load gatk/4.0-8-Java-1.8.0_92
      gatk -Xmx32G ApplyVQSR -O tmp.indel.recalibrated.vcf -V ${vcf} --recal-file ${indel_recal} --tranches-file ${indel_tranches} \
      --truth-sensitivity-filter-level ${params.indel_filter_level} -create-output-variant-index true -mode INDEL

      gatk -Xmx32G ApplyVQSR -O ${recal_vcf} -V tmp.indel.recalibrated.vcf --recal-file ${snv_recal} --tranches-file ${snv_tranches} \
      --truth-sensitivity-filter-level ${params.snp_filter_level} --create-output-variant-index true -mode SNP
    """
}


process runCollectVariantCallMetrics {
  tag "$BATCH"
  label "mediumcpu"
  conda "picard"
  publishDir "${OUTDIR}/${BATCH}/VariantCallingStats"
  input:
    tuple file(recal_vcf), file(recal_vcf_index) from runApplyRecalibrationOut
  output:
    file("${BATCH}.variant_calling_detail_metrics") into runCollectVariantCallMetricsOut
  script:
    """
    picard -Xmx24g CollectVariantCallingMetrics INPUT=${recal_vcf} OUTPUT=${BATCH} DBSNP=${DBSNP} 
    """
}

process runNormalization {
  tag "$BATCH"
  label "mediumcpu"
  conda "bcftools"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  input:
    tuple file(recal_vcf), file(recal_vcf_index) from runNormalizationIn 
  output:
    tuple file(outVcf), file(outVcf_index) into runNormalizationOut
  script:
    outVcf = BATCH + ".norm.vcf.gz"
    outVcf_index = BATCH + ".norm.vcf.gz.tbi"
    """
      bcftools norm -m -any --threads 4 -Oz -o tmp.norm.vcf.gz ${recal_vcf}
      tabix -p vcf tmp.norm.vcf.gz
      bcftools view -e "AC=0" -Oz -o ${outVcf} tmp.norm.vcf.gz
      tabix -p vcf ${outVcf}
    """
}


process runSplitByChr {
  tag "$BATCH|$chr"
  publishDir "${OUTDIR}/${BATCH}/VCF/byChr"
  conda "bcftools"
  input:
    tuple file(vcf), file(vcf_index) from runNormalizationOut
    each chr from chr_vep
  output:
    tuple chr, file(outVcf), file(outVcf_index) into runSplitByChrOut
  script:
    outVcf = BATCH + ".chr" + chr + ".norm.vcf.gz"
    outVcf_index = BATCH + ".chr" + chr + ".norm.vcf.gz.tbi"
    """
      bcftools view -r ${chr} -Oz -o ${outVcf} ${vcf}
      tabix -p vcf ${outVcf}
    """
}

process runVEP {
  tag "$BATCH|$chr"
  label "highcpu"
  publishDir "${OUTDIR}/${BATCH}/VCF/byChr"
  input:
    tuple chr, file(vcf), file(vcf_index) from runSplitByChrOut
  output:
    file(outVcf) into runVEPVcfOut
    file(outVcf_index) into runVEPVcfIndexOut
  script:
    outVcf = BATCH + ".chr" + chr + ".vep.vcf.gz"
    outVcf_index = BATCH + ".chr" + chr + ".vep.vcf.gz.tbi"
    """
      module load VEP/89-goolf-1.7.20
      \${EBROOTVEP}/vep --fork 8 --dir \${EBROOTVEP}/.vep  --assembly GRCh38 \
      --everything --vcf --allele_number --no_stats --cache --offline --force_overwrite --compress_output bgzip \
      --plugin dbNSFP,${dbNSFP_dir}/dbNSFP.gz,Polyphen2_HVAR_pred,CADD_phred,SIFT_pred,FATHMM_pred,MutationTaster_pred,MetaSVM_score,MetaSVM_pred,MetaLR_score,MetaLR_pred,Reliability_index \
      --plugin LoF,human_ancestor_fa:\${EBROOTVEP}/.vep/Plugins/loftee/human_ancestor.fa.gz,filter_position:0.05,min_intron_size:15 \
      -i ${vcf} -o ${outVcf}
      tabix -p vcf ${outVcf}
    """
}



process runCombineVEPVcf {
  tag "$BATCH"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  label "mediumcpu"
  conda "gatk4"
  input:
    file(vcf) from runVEPVcfOut.collect()
    file(vcf_index) from runVEPVcfIndexOut.collect()
  output:
    tuple file(outVcf), file(outVcf_index) into runCombineVEPVcfOut
  script:
    outVcf = BATCH + ".vep.vcf.gz"
    outVcf_index = BATCH + ".vep.vcf.gz.tbi"
    """
      gatk -Xmx16G SortVcf --INPUT ${vcf.join(" --INPUT ")} --OUTPUT ${outVcf} --SEQUENCE_DICTIONARY ${REF_DIC}
    """
}


workflow.onComplete {
  log.info "========================================="
  log.info "Duration:           $workflow.duration"
  log.info "========================================="

  def msg = """\
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    .stripIndent()

  sendMail(to: 'huangk3@nih.gov', subject: 'WES pipeline execution', body: msg)
}

// ------------------------------------------------------------------------------------------------------------
//
// Read input file and save it into list of lists
//
// ------------------------------------------------------------------------------------------------------------
def logParams(p, n) {
  File file = new File(n)
  file.write "Parameter:\tValue\n"

  for(s in p) {
     file << "${s.key}:\t${s.value}\n"
  }
}

