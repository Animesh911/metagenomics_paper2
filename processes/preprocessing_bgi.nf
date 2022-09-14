process First_Fastqcb {
    publishDir path: { "$baseDir/preprocessing/bgi/${sample_id}/fastqc" }
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


process Second_Fastqcb {
  publishDir "$baseDir/preprocessing/bgi/${sample_id}/fastqc"
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

process Multiqcb {
  publishDir "$baseDir/preprocessing/bgi/"
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


process Clumpifyb {
  publishDir "$baseDir/preprocessing/bgi/${sample_id}/clumpify"
  tag "$sample_id"

  input:
  tuple val(sample_id), path (reads)
  path r3
  
  output:
  tuple val(sample_id), path ("*.fq.gz"), emit: trim
  
  shell:
  """
  bbduk.sh -Xmx20g in1=${reads[0]} in2=${reads[1]} out1=${reads[0].getSimpleName()}.fq.gz out2=${reads[1].getSimpleName()}.fq.gz ref=!{r3} stats=Stats_${reads[0].getSimpleName()}.txt minlen=51 qtrim=r trimq=10 forcetrimleft=3
  """
}

process Repairb {
  publishDir "$baseDir/preprocessing/bgi/${sample_id}/repair"
  container 'hunter:latest'
  tag "$sample_id"
  
  input:
  tuple val(sample_id), path (r1)

  output:
  tuple val(sample_id), path ("*_fixed.fq.gz"), emit: synced_r1

  shell:
  """
  repair.sh -Xmx20g in=!{r1[0]} in2=!{r1[1]} out=!{r1[0].getSimpleName()}_fixed.fq.gz out2=!{r1[1].getSimpleName()}_fixed.fq.gz outs=!{r1[0].getSimpleName()}_singletons.fq.gz repair
  """
}
