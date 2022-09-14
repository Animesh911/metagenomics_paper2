process First_Fastqc {
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/fastqc"
  tag "Quality checking raw reads: $sample_id"

  input:
  tuple val(sample_id), path(reads)

  output:
  path "${sample_id}_raw"

  script:
  """
  mkdir ${sample_id}_raw
  fastqc -o ${sample_id}_raw -f fastq -q ${reads}
  """
}

process Second_Fastqc {
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/fastqc"
  tag "Quality checking preprocessed reads: $sample_id"
  
  input:
  tuple val(sample_id), path (reads)

  output:
  path "${sample_id}_post"

  script:
  """
  mkdir ${sample_id}_post
  fastqc -o ${sample_id}_post -f fastq -q ${reads}
  """
}

process Multiqc {
  publishDir "$baseDir/preprocessing/illumina/"
  tag "Generating multiQC report"

  input:
  file(fastqc_out)
  file(fastqc_out)

  output:
  path ("multiqc"), emit: dir

  shell:
  '''
  multiqc --outdir multiqc .
  '''
}

process Optical {
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/optical"
  tag "$sample_id"

  input:
  tuple val(sample_id), path(reads)
  
  output:
  tuple val(sample_id), path ("*_clump.fastq.gz"), emit: optical
  
  shell:
  """
  clumpify.sh -Xmx20g in1=${reads[0]} in2=${reads[1]} out1=${reads[0].getSimpleName()}_clump.fastq.gz out2=${reads[1].getSimpleName()}_clump.fastq.gz dedupe optical dist=40
  """
}

process Clumpify {
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/clumpify"
  tag "$sample_id"

  input:
  tuple val(sample_id), path (reads)
  path r3
  
  output:
  tuple val(sample_id), path ("*_bbduk.fastq.gz"), emit: trim_r1
  
  shell:
  """
  bbduk.sh -Xmx20g in1=${reads[0]} in2=${reads[1]} out1=${reads[0].getSimpleName()}_bbduk.fastq.gz  out2=${reads[1].getSimpleName()}_bbduk.fastq.gz ref=!{r3} forcetrimleft=17 ktrim=r minlen=51 qtrim=r trimq=10 tbo=t mink=11 hdist=1
  """
}

process Fastqscreen {
  container 'hunter:latest'
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/fastq_screen"
  tag "$sample_id"

  input:
  tuple val(sample_id), path (r1)
  path r3

  output:
  tuple val(sample_id), path ("*.tagged_filter.fastq.gz"), emit: contamin_r1

  shell:
  """
  fastq_screen --nohits !{r1[0]} !{r1[1]} --conf !{r3} --aligner bowtie2
  """
}

process Repair {
  publishDir "$baseDir/preprocessing/illumina/${sample_id}/repair"
  container 'hunter:latest'
  tag "$sample_id"
  
  input:
  tuple val(sample_id), path (r1)

  output:
  tuple val(sample_id), path ("*.desensitized.fastq.gz"), emit: synced_r1

  shell:
  """
  repair.sh -Xmx20g in=!{r1[0]} in2=!{r1[1]} out1=!{r1[0].getSimpleName()}.desensitized.fastq.gz  out2=!{r1[1].getSimpleName()}.desensitized.fastq.gz
  """
}

