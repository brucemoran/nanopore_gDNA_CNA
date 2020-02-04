#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '--------------------------------------------------------------'
  log.info 'NEXTFLOW NANOPORE FASTQ, ALIGN, QC AND CNA ANALYSIS WITH GUPPY'
  log.info '--------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run nanopore_gDNA_CNA.nf \
              --dataDir "path/to/data/20181010_1732_MMUH_Nanopore_101018" \
              --barcodeIDmap "barcode.sampleID.map.csv" \
              --fa hg38.fa'

  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --dataDir     STRING      path from {workflow.launchDir} where Minion data found in fast5 format; NB can have multiple sub-directories (0,1,2 etc) containing fast5s'
  log.info '    --barcodeIDMap      STRING      csv file of barcode,sampleID (with that as header)'
  log.info '    --fa      STRING      path to fasta'
  log.info ''
  exit 1
}

/* 0.0: Global Variables
*/

/* General
*/
params.outDir = "${workflow.launchDir}/analysis"
params.dataDir = "${workflow.launchDir}/data"

/* 0.0: minimap2 index
*/
process index {

  output:
  file('minimap.mmi') into mnm_index
  file('fasta.fa') into pbsv_fasta
//  file("fasta.fa.gz") into medaka_fasta

  script:
  """
  {
  if [[ ${params.fa} =~ ".gz\$" ]];then
    cp ${params.fa} fasta.fa.gz
    gunzip -c fasta.fa.gz > fasta.fa
    gunzip -c fasta.fa.gz | bgzip -c > fasta.fa.bgz
    samtools faidx fasta.fa.bgz
  else
    cp ${params.fa} fasta.fa
    gzip -c fasta.fa > fasta.fa.gz
    gunzip -c fasta.fa.gz | bgzip -c > fasta.fa.bgz
    samtools faidx fasta.fa.bgz
  fi

  minimap2 -x map-ont -d minimap.mmi fasta.fa.gz
  } 2>&1 | tee index.log.txt
  """
}

/* 1.0: QC sequencing run
*
*/
process guppy {

  publishDir "$params.outDir/guppy", mode: "copy"

  output:
//  file('guppy/*fastq') into fastq
//  file('guppy') into medaka_fastq
  set file('guppy/sequencing_summary.txt'), file('guppy/barcoder/barcoding_summary.txt') into qc_summary
  file('guppy/barcoder') into barcodedir

  script:
  """
  TMEM=\$(echo ${task.memory} | sed 's/ GB//')
  ##callers is cpu; threads-per-caller is mem/cpu
  TPCS=\$(( \$TMEM / ${task.cpus} ))
  guppy_basecaller \
    -s guppy \
    -c ${params.guppyConfig} \
    -i ${params.dataDir} \
    -r \
    --cpu_threads_per_caller \$TPCS \
    --num_callers ${task.cpus} \
    --min_qscore 7 \
    -q 0

  guppy_barcoder \
    -s guppy/barcoder \
    -i guppy \
    --worker_threads ${task.cpus} \
    -q 0
  """
}

// process guppy_f5 {
//
//   publishDir "$params.outDir/guppy", mode: "copy"
//
//   output:
//   file('guppy') into medaka_fast5
//
//   script:
//   """
//   ##work out how many callers and how much threads
//   TMEM=\$(echo ${task.memory} | sed 's/ GB//')
//   ##callers is cpu; threads-per-caller is mem/cpu
//   TPCS=\$(( \$TMEM / ${task.cpus} ))
//
//   guppy_basecaller \
//     --fast5_out \
//     --post_out \
//     -s guppy \
//     -c ${params.guppyConfig} \
//     -i ${params.dataDir} \
//     -r \
//     --cpu_threads_per_caller \$TPCS \
//     --num_callers ${task.cpus} \
//     --min_qscore 7 \
//     --flowcell ${params.flowcell} \
//     --kit ${params.kit} \
//     -q 0
//   """
// }

/* 1.01 Medaka Methylation Calling
*/
// process medaka {
//
//   publishDir "$params.outDir/medaka", mode: "copy"
//
//   input:
//   file(fa) from medaka_fasta
//   file(f5) from medaka_fast5
//
//   output:
//   file('*') into medout
//
//   script:
//   """
//   {
//     OUTBAM=meth.bam
//     medaka methylation guppy2sam $f5 $fa \
//         --workers ${task.cpus} --recursive | \
//     samtools sort -@ ${task.cpus} | \
//     samtools view -b -@ ${task.cpus} > medaka.meth.bam
//   } 2>&1 | tee > "medaka.log.txt"
//   """
// }

/* 1.1: QC sequencing run
*
*/
process pycoqc {

  publishDir "$params.outDir/pycoQC", mode: "copy"

  input:
  set file(seqsummary), file(barcodesummary) from qc_summary

  output:
  file('*') into completed_11

  script:
  """
  pycoQC \
    --file $seqsummary \
    --barcode_file $barcodesummary \
    --outfile pycoQC.html
  """
}

/* 2.0: Alignment
* take in $params.barcodeIDMap, held in $params.dataDir
*/
Channel.fromPath("$params.barcodeIDmap", type: 'file')
       .splitCsv( header: true )
       .set { barcodeIDmap }

process minimap2 {

  publishDir "$params.outDir/$sampleID/minimap2", mode: 'copy', pattern: '*[!.gz,.txt]'
  publishDir "$params.outDir/$sampleID/fastq", mode: 'copy', pattern: '*.fastq.gz'

  input:
  set val(barcode), val(sampleID) from barcodeIDmap
  each file(barcodebases) from barcodedir
  each file(mmi) from mnm_index

  output:
  set val(sampleID), file('*.bam') into (nanostat_bam, pbsv_bam)
  file('*.guppy.pass.fastq.gz') into mnm_fastq

  script:
  """
  INPUT=\$(find $barcodebases/$barcode | grep fastq)
  cp \$INPUT "${sampleID}.guppy.pass.fastq"
  gzip "${sampleID}.guppy.pass.fastq"
  DAT=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ONT\\tSM:$sampleID\\tDS:ONT\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"
  minimap2 \
    -t ${task.cpus} \
    -ax map-ont \
    -R \$RGLINE \
    -Y \
    $mmi \
    "${sampleID}.guppy.pass.fastq.gz" | \
  samtools sort -T "tmp."$sampleID - | \
  samtools view -Shb - > $sampleID".bam"
  """
}

process nanostat {

  publishDir "$params.outDir/$sampleID/nanostat", mode: 'copy', pattern: '*.nanostat.txt'

  input:
  set val(sampleID), file(bam) from nanostat_bam

  output:
  file('*.nanostat.txt') into nnosttc

  script:
  """
  NanoStat --bam $bam > $sampleID".nanostat.txt"
  """
}

process pbsv_discover_call {

  publishDir "$params.outDir/$sampleID/pbsv", mode: "copy"

  input:
  file(fasta) from pbsv_fasta
  set val(sampleID), file(bam) from pbsv_bam

  output:
  set val(sampleID), file("${sampleID}.svsig.vcf") into pbsv_vcf
  file('*') into pbsv_complete

  script:
  """
  pbsv discover $bam ${sampleID}.svsig.gz
  pbsv call $fasta ${sampleID}.svsig.gz ${sampleID}.svsig.vcf
  """
}

process pbsv_vcf {

  publishDir "$params.outDir/$sampleID/pbsv", mode: "copy"

  input:
  set val(sampleID), file("${sampleID}.svsig.vcf") from pbsv_vcf

  output:
  file('*') into pbsv_complete2

  script:
  """
  pbsv discover $bam ${sampleID}.svsig.gz
  pbsv call $fasta ${sampleID}.svsig.gz ${sampleID}.svsig.vcf
  """
}

// process fctcsv {
//
//   label 'c8_24G_cpu_mem'
//
//   publishDir "$params.outDir/$caseID/$sampleID/facets"
//
//   input:
//   set val(sampleID), file(bam), file(bai), val(germlineID), file(germlinebam), file(germlinebai) from facetsing
//   set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
//   file(facetsR) from facetscallscript
//
//   output:
//   val(sampleID) into completed2_1
//   file('*') into facetsoutputR
//
//   script:
//   """
//   CSVFILE=\$(echo $bam | sed 's/bam/facets.r10.csv/')
//
//   {
//     snp-pileup \
//       $dbsnp \
//       -r 10 \
//       -p \
//       \$CSVFILE \
//       $germlinebam \
//       $bam
//
//     Rscript --vanilla $facetsR \$CSVFILE
//   } 2>&1 | tee > $sampleID".facets_snpp_call.log.txt"
//   """
// }
