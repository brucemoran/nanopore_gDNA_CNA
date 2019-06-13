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
              --dataDir "/full/path/to/data/ 20181010_1732_MMUH_Nanopore_101018" \
              --barcodeIDmap "barcode.sampleID.map.csv" \
              -c "nanopore_gDNA_CNA.nextflow.config" \
              -with-report "nanopore_gDNA_CNA.report.html" \
              -with-timeline "nanopore_gDNA_CNA.timeline.html"'

  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --dataDir     STRING      relative path from baseDir where Minion data found in fast5 format; NB can have multiple sub-directories (0,1,2 etc) containing fast5s'
  log.info '    --barcodeIDMap      STRING      csv file of barcode,sampleID (with that as header)'
  log.info ''
  exit 1
}

/* 0.0: Global Variables
*/

/* General
*/
params.outDir = "$baseDir/analysis"
params.dataDir = "$baseDir/data"

/* 0.0: minimap2 index
*/
if (!params.fa){
  params.fain = params.hg38fa
}
else {
  params.fain = params.fa
}


process index {

  label 'c5_15G_cpu_mem'

  output:
  file('minimap.mmi') into mnm_index
  file('*.fa') into cnvnator_fasta

  script:
  """
  wget ${params.fain}
  FNA=\$(basename ${params.fain})
  minimap2 -x map-ont -d minimap.mmi \$FNA

  ##for CNVnator
  CHRS=\$(echo chr{1..21})
  for CHR in chr{1..21};
    do echo \$CHR > 1;
    samtools faidx -r 1 \$FNA > \$CHR".fa";
  done
  """
}

/* 1.0: QC sequencing run
*
*/
process guppy {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/guppy", mode: "copy"

  output:
  file('guppy/*fastq') into fastq
  set file('guppy/sequencing_summary.txt'), file('guppy/barcoder/barcoding_summary.txt') into qc_summary
  file('guppy/barcoder') into barcodedir

  script:
  """
  guppy_basecaller \
    -s guppy \
    -c dna_r9.4.1_450bps_flipflop.cfg \
    -i ${params.dataDir} \
    -r \
    --cpu_threads_per_caller ${task.cpus} \
    --num_callers 1 \
    --min_qscore 7 \
    # --flowcell ${params.flowcell} \
    # --kit ${params.kit} \
    -q 0

  guppy_barcoder \
    -s guppy/barcoder \
    -i guppy \
    --worker_threads ${task.cpus} \
    -q 0
  """
}

/* 1.1: QC sequencing run
*
*/
process pycoqc {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/pycoQC", mode: "copy"

  input:
  set file(seqsummary), file(barcodesummary) from qc_summary

  output:
  file('*') into completed_11

  script:
  """
  pycoQC \
    --summary_file $seqsummary \
    --barcode_file $barcodesummary \
    --html_outfile pycoQC.html
  """
}

/* 2.0: Alignment
* take in $params.barcodeIDMap, held in $params.dataDir
*/
Channel.fromPath("$params.barcodeIDmap", type: 'file')
       .splitCsv( header: true )
       .set { barcodeIDmap }

process minimap2 {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/minimap2", mode: 'copy', pattern: '*[!.gz,.txt]'
  publishDir "$params.outDir/$sampleID/fastq", mode: 'copy', pattern: '*.fastq.gz'

  input:
  set val(barcode), val(sampleID) from barcodeIDmap
  each file(barcodebases) from barcodedir
  each file(mmi) from mnm_index

  output:
  set val(sampleID), file('*.bam') into (nanostat, cnvnator)
  file('*.guppy.pass.fastq.gz') into mnm_fastq

  script:
  """
  INPUT=\$(find $barcodebases/$barcode -name *fastq)
  cp \$INPUT "${sampleID}.guppy.pass.fastq"
  gzip "${sampleID}.guppy.pass.fastq"

  minimap2 \
    -t ${task.cpus} \
    -ax map-ont \
    $mmi \
    "${sampleID}.guppy.pass.fastq.gz" | \
  samtools sort -T "tmp."$sampleID - | \
  samtools view -Shb - > $sampleID".bam"
  """
}

process nanostat {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$sampleID/nanostat", mode: 'copy', pattern: '*.nanostat.txt'

  input:
  set val(sampleID), file(bam) from nanostat

  output:
  file('*.nanostat.txt') into nnosttc

  script:
  """
  NanoStat --bam $bam > $sampleID".nanostat.txt"
  """
}

BINS = Channel.from(5000, 10000, 50000, 100000)

process cnvnator {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$sampleID/CNVnator", mode: "copy"

  input:
  each binsize from BINS
  each fasta from cnvnator_fasta
  set val(sampleID), file(bam) from cnvnator

  output:
  file('*') into cnvnator_complete

  script:
  """

  cnvnator \
    -root $sampleID".root" \
    -chrom \$CHRS \
    -tree $bam

  cnvnator \
    -root $sampleID".root" \
    -chrom \$CHRS \
    -his $binsize \
    -d ./

  cnvnator \
    -root $sampleID".root" \
    -chrom \$CHRS \
    -stat $binsize

  cnvnator \
    -root $sampleID".root" \
    -chrom \$CHRS \
    -partition $binsize

  cnvnator \
    -root $sampleID".root" \
    -chrom \$CHRS \
    -call $binsize > $sampleID".cnvnator.out.txt"

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
