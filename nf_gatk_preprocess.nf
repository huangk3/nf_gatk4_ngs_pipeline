#!/bin/env nextflow

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

logParams(params, "nextflow_parameters.txt")

VERSION = "0.1"

// Header log info
log.info "========================================="
log.info "GRIS WES pipeline for Exome-Seq Preprocessing v${VERSION}"
log.info "Nextflow Version:     $workflow.nextflow.version"
log.info "Command Line:         $workflow.commandLine"
log.info "========================================="

Channel.fromFilePairs("/nethome/huangk3/bcbb/Andrew_Demo/WESTrainingCourse/daughter*/*{1,2}.fq").map{item -> item[1]}.into{fastqPairs; extractSampleSizes}

N=extractSampleSizes.count()

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
  publishDir "${OUTDIR}/${BATCH}/${sampleID}/FastqToSam/"
  input: 
    tuple file(fq1), file(fq2) from fastqPairs

  output:
    tuple sampleID, file(outBam), file(outBai) into runFastqToBamOut
    val sampleID into sampleID_list, sampleID_mapping
  script:

  sampleID = fq1.getSimpleName()
  outBam = sampleID + ".sorted.bam"
  outBai = sampleID + ".sorted.bai"
  """
    bwa mem -t 8 -R "@RG\\tID:${sampleID}\\tLB:${sampleID}\\tSM:${sampleID}\\tPL:ILLUMINA" ${REF}  ${fq1} ${fq2} | \
    java -Xmx30g -jar \${EBROOTPICARD}/picard.jar SortSam I=/dev/stdin O=${outBam} \
    MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate
  """

}


process runMarkDuplicates {
  tag "${BATCH}|${sampleID}"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/${sampleID}/MarkDuplicate/"
  input:
    tuple sampleID, file(bam), file(bai) from runFastqToBamOut

  output:
    tuple sampleID, file(outBam) , file(outBai), file(outMetrics) into runMarkDuplicateOut
    tuple BATCH, file(outBam) into runFastQC
  script:
    outBam = sampleID + ".markdup.bam"
    outBai = sampleID + ".markdup.bai"
    outMetrics = sampleID + ".markdup.metrics.txt" 
  """
    java -jar \${EBROOTPICARD}/picard.jar MarkDuplicates I=${bam} O=${outBam} \
    VALIDATION_STRINGENCY=LENIENT AS=true CREATE_INDEX=true M=${outMetrics} 
  """
}


process runBaseRecalibrator {
  tag "$BATCH|$sampleID"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/${sampleID}/BQSR/"
  
  input:
    tuple sampleID, file(bam) , file(bai), file(mkdup_metrics) from runMarkDuplicateOut

  output:
    tuple sampleID, file(bam) , file(bai), file(recal_table) into runBaseRecalibratorOut

  script:
    recal_table=sampleID + ".recal_data.table"
    
    """
      \${EBROOTGATK}/gatk BaseRecalibrator -R ${REF} -I ${bam} -O ${recal_table} --known-sites ${DBSNP} --known-sites ${OneThousand} -known-sites ${MILLS}
    """
}

process runApplyBQSR {
  tag "$BATCH|$sampleID"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/${sampleID}/BQSR/"

  input:
    tuple sampleID, file(bam) , file(bai), file(recal_table) from runBaseRecalibratorOut

  output:
    tuple sampleID, file(outBam), file(outBai) into runApplyBQSROut
   
  script:
    outBam = sampleID + ".bqsr.bam"
    outBai = sampleID + ".bqsr.bai"

    """
      \${EBROOTGATK}/gatk ApplyBQSR -R ${REF} -I ${bam} -bqsr ${recal_table} -O ${outBam}
    """
}

process runHaplotypeCaller {
  tag "$BATCH|$sampleID"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/${sampleID}/HaplotypeCaller/"

  input:
    tuple sampleID, file(bam) , file(bai) from runApplyBQSROut

  output:
    tuple BATCH, file(gvcf), file(gvcf_index) into runHaplotypeCallerOut

  script:
    gvcf = sampleID + ".g.vcf.gz"
    gvcf_index = sampleID + ".g.vcf.gz.tbi"

    """
      \${EBROOTGATK}/gatk HaplotypeCaller -R ${REF} -I ${bam} -O ${gvcf} -stand-call-conf 20.0 --dbsnp ${DBSNP} -L ${TARGETS} -ERC GVCF
      tabix -p vcf -f ${gvcf}
    """
}

sampleID_list.collectFile(newLine: true, storeDir: "${OUTDIR}/${BATCH}") {item->
  ["${BATCH}.sampleList.txt", item]
}

sampleID_mapping.collectFile(newLine: true, storeDir: "${OUTDIR}/${BATCH}") {item->
  ["${BATCH}.gVCF.mapping.txt", item + "\t" + "${OUTDIR}/${BATCH}/"+ item +"/HaplotypeCaller/" + item + ".g.vcf.gz"]
}

//runHaplotypeCallerOut_grouped_by_batch = runHaplotypeCallerOut.groupTuple(by: [0])
total_intervals = Channel.fromPath("${TARGETS}").countLines()
process runDynamicallyCombineIntervals {
  publishDir "${OUTDIR}/${BATCH}/splittedIntervals/"
  input:
    val n from N
    val total from total_intervals

  output:
    file("p*.intervals") into runDynamicallyCombineIntervalsOut

  script:
    lines_per_file =  (1 + total / (n + 30)).toInteger()
    """
      split -d -a4 -l $lines_per_file --additional-suffix ".intervals" ${TARGETS} "p"
    """
}


process runGenomicsDBImport {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/GenomicsDB"

  input:
    tuple BATCH, file(gvcf), file(gvcf_index) from runHaplotypeCallerOut
    each file(interval) from runDynamicallyCombineIntervalsOut
  output:
    tuple BATCH, p, file(interval), file(genomeDb) into runGenomicsDBImportOut

  script:
    p = interval.getSimpleName()
    genomeDb = BATCH + "." + p + ".genomeDb"
  """
    \$EBROOTGATK/gatk GenomicsDBImport --genomicsdb-workspace-path ${genomeDb} \
    --overwrite-existing-genomicsdb-workspace true --batch-size 50 -ip ${params.ip} -L ${interval} --sample-name-map ${OUTDIR}/${BATCH}/${BATCH}.gVCF.mapping.txt
  """ 
}

process runGenotypeGVCFs {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/vcfByIntervals"

  input:
    tuple BATCH, p, file(interval), file(genomeDb) from runGenomicsDBImportOut

  output:
    tuple BATCH, p, file(outVcf), file(outVcf_index) into runGenotypeGVCFsOut

  script:
    outVcf = BATCH + "." + p + ".vcf.gz"
    outVcf_index = BATCH + "." + p + ".vcf.gz.tbi"
  """
    \$EBROOTGATK/gatk GenotypeGVCFs -R ${REF} -V gendb://${genomeDb} -O ${outVcf} -G StandardAnnotation --only-output-calls-starting-in-intervals \
    --use-new-qual-calculator -L ${interval}
  """
}

runGenotypeGVCFsOut_group_by_batch=runGenotypeGVCFsOut.groupTuple(by: [0])

process runSortVcfs {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/VCF"

  input:
    tuple BATCH, p_list, vcf_list, vcf_index_list from runGenotypeGVCFsOut_group_by_batch
  output:
    tuple file(outVcf), file(outVcf_index) into runSortVcfsOut, ApplyVQSR  
  script:
    outVcf = BATCH + ".vcf.gz"
    outVcf_index = BATCH + ".vcf.gz.tbi"
    
  """
    \$EBROOTGATK/gatk SortVcf --INPUT ${vcf_list.join(" --INPUT ")} --OUTPUT ${outVcf} --SEQUENCE_DICTIONARY ${REF_DIC}
  """
} 

process runHardFilterAndMakeSitesOnlyVcf {
  tag "$BATCH"
  label "mediumcpu"
  publishDir "${OUTDIR}/${BATCH}/VCF"
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
      \$EBROOTGATK/gatk VariantFiltration --filter-expression "ExcessHet > ${params.excess_het_threshold}" --filter-name ExcessHet -O ${filtered_vcf} -V ${vcf}
      java -jar \${EBROOTPICARD}/picard.jar MakeSitesOnlyVcf INPUT= ${filtered_vcf} OUTPUT=${sitesOnly_vcf}
    """
}

process runIndelRecal {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/VQSR"
  input:
    tuple file(vcf), file(vcf_index) from indelRecal
  output:
    tuple file(indel_recal), file(indel_recal_index), file(indel_tranches) into runIndelRecalOut
  script:
    indel_recal = BATCH + ".indel.recal"
    indel_recal_index = BATCH + ".indel.recal.idx"
    indel_tranches = BATCH + ".indel.tranches"
  """
    \$EBROOTGATK/gatk VariantRecalibrator -R ${REF} -V ${vcf} --output ${indel_recal} --tranches-file ${indel_tranches} --trust-all-polymorphic \
    -an ${indel_recal_anno_values.join(" -an ")} -tranche ${indel_recal_tranche_values.join(" -tranche ")} \
    -mode INDEL --max-gaussians 4 -resource mills,known=false,training=true,truth=true,prior=12:${MILLS} \
    -resource axiomPoly,known=false,training=true,truth=false,prior=10:${AXIOM} \
    -resource dbsnp,known=true,training=false,truth=false,prior=2:${DBSNP}
  """
}

process runSnvRecal {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/VQSR"
  echo true
  input:
    tuple file(vcf), file(vcf_index) from snvRecal
  output:
    tuple file(snv_recal), file(snv_recal_index), file(snv_tranches) into runSnvRecalOut
  script:
    snp_recal = BATCH + ".snv.recal"
    snp_recal_index = BATCH + ".snv.recal.idx"
    snp_tranches = BATCH + ".snv.tranches"
  """
    \$EBROOTGATK/gatk VariantRecalibrator -R ${REF} -V ${vcf} --output ${snp_recal} --tranches-file ${snp_tranches} --trust-all-polymorphic \
    -an ${snp_recal_anno_values.join(" -an ")} -tranche ${snp_recal_tranche_values.join(" -tranche ")} \
    -mode SNP --max-gaussians 6 \
    -resource hapmap,known=false,training=true,truth=true,prior=15:${HAPMAP} \
    -resource omni,known=false,training=true,truth=true,prior=12:${OMINI} \
    -resource 1000G,known=false,training=true,truth=false,prior=10:${OneThousand} \
    -resource dbsnp,known=true,training=false,truth=false,prior=7:${DBSNP}
  """
}


process runApplyVQSR {
  tag "$BATCH"
  label "highmem"
  publishDir "${OUTDIR}/${BATCH}/VCF"
  input:
    tuple file(indel_recal), file(indel_recal_index), file(indel_tranches) from runIndelRecalOut
    tuple file(snv_recal), file(snv_recal_index), file(snv_tranches) from runSnvRecalOut
    tuple file(vcf), file(vcf_index) from ApplyVQSR
  output:
    tuple file(recal_vcf), file(recal_vcf_index) into runApplyRecalibrationOut
  script:
    recal_vcf = BATCH + ".recal.vcf.gz"
    recal_vcf_index = BATCH + ".recal.vcf.gz.tbi"
    """
      \$EBROOTGATK/gatk ApplyVQSR -O tmp.indel.recalibrated.vcf -V ${vcf} --recal-file ${indel_recal} --tranches-file ${indel_tranches} \
      --truth-sensitivity-filter-level ${params.indel_filter_level} -create-output-variant-index true -mode INDEL

      \$EBROOTGATK/gatk ApplyVQSR -O ${recal_vcf} -V tmp.indel.recalibrated.vcf --recal-file ${snv_recal} --tranches-file ${snv_tranches} \
      --truth-sensitivity-filter-level ${params.snp_filter_level} --create-output-variant-index true -mode SNP
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

