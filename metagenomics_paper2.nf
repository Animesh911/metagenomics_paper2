#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
===================================================================================
                    Marine Metagenomics analysis Pipeline                     
===================================================================================

nextflow/metagenomics Analysis Pipeline Paper2. Started 2022-09-12
#### Homepage / Documentation
https://github.com/Animesh911/Metagenomics_paper2
#### Author
Animesh Kumar <animesh.kumar@uit.no>
TOKEN: eyJ0aWQiOiA1MjcxfS42YTExNWVjNzFiYWM1MjliYzFmY2YxNTNhNzE4MzgxZTA1Y2FjOGIy
-----------------------------------------------------------------------------------
*/

// Help Message
def helpMessage() {
    log.info """
    *********Marine metagenomics analysis*********
    
    Usage:
    nextflow run metagenomics_paper2.nf -with-docker --illumina /home/aku048/storage/paper2/raw_reads1/illumina --bgi /home/aku048/storage/paper2/raw_reads1/bgi --outputill /home/aku048/storage/paper2/preprocessing/illumina/ --outputbgi /home/aku048/storage/paper2/preprocessing/bgi/ -with-tower -resume
    
    Input:
    --bgi                      path to the directory containing the BGI read file (fastq.gz) (default: $params.bgi)
    --illumina                 path to the directory containing the illumina read file (fastq.gz) (default: $params.illumina)
    
    Output:
    --outputbgi                path to the output directory (default: $params.outputbgi)
    --outputill                path to the output directory (default: $params.outputill)
    
    Taxonomic database:
    --kaiju_refseq             Database (Refseq bactarch) for taxonomic binning with kaiju (default: $params.refseqkaiju). E.g. "/storage/aku048/database/refseq/kaiju1.7.3/refseq_arch_bact_default/refseq/kaiju_db_refseq.fmi"
    --kaiju_refseq_nodes       Database (Refseq bactarch) for taxonomic nodes with kaiju bactarch (default: $params.refseqkaijunodes). E.g. "/storage/aku048/database/refseq/kaiju1.7.3/refseq_arch_bact_default/nodes.dmp"
    --kraken_refseq            Database (Refseq bactarch) for taxonomic binning with kraken2 (default: $params.refseqkraken). E.g. "/storage/aku048/database/refseq/kraken2/"
    
    --kaiju_marRef             Database (marRef prot) for taxonomic binning with kaiju (default: $params.mar_refprot). E.g. "/storage/aku048/database/mar5/kaiju_mar_ref_protein/marref_proteins_V5.fmi"
    --kraken_marRef            Database (marRef nucl) for taxonomic binning with kraken2 (default: $params.mar_refnucl). E.g. "/storage/aku048/database/mar5/kraken_mar_ref_nucl/database/"
    --kaiju_marDb              Database (marDb prot) for taxonomic binning with kaiju (default: $params.mar_dbprot). E.g. "/storage/aku048/database/mar5/kaiju_mar_db_protein/masked_segmar5_db_prot.fmi"
    --kraken_marDb             Database (marDb nucl) for taxonomic binning with kraken2 (default: $params.mar_dbnucl). E.g. "/storage/aku048/database/mar5/kraken_mar_db_nucl/database/"
    
    --kaiju_marRefDb           Database (marRefDb combined prot) for taxonomic binning with kaiju (default: $params.mar_ref_dbprot). E.g. "/storage/aku048/database/mar5/kaiju_mar_ref_db_protein/mar5_ref_db_prot.fmi"
    --kraken_marRefDb          Database (marRefDb combined nucl) for taxonomic binning with kraken2 (default: $params.mar_ref_dbnucl). E.g. "/storage/aku048/database/mar5/kraken_mar_ref_db_nucl/database/"
    --help                     This usage statement (default: $params.help=false)
    
"""
}

// Show help message
if (params.help) {
	helpMessage()
	exit 0
}


//params.reads = "$baseDir/raw_data/ST192_{fw,rv}.fastq.gz"
params.adapter_ill = "$baseDir/adapters/adapters.fa"
params.adapter_bgi = "$baseDir/adapters/BGI.fa"
params.multiqc = "$baseDir/multiqc"
params.fastqscreen = "$baseDir/fastq_screen.conf"
params.help = false

/*
//Default Database
params.refseqkaiju = "/storage/aku048/database/refseq/kaiju1.7.3/refseq/kaiju_refseq_bact_arch.fmi"
//params.refseqkaiju = "/storage/aku048/database/refseq/kaiju1.7.3/refseq_arch_bact_default/refseq/kaiju_db_refseq.fmi"
params.refseqkaijunodes = "/storage/aku048/database/refseq/kaiju1.7.3/refseq_arch_bact_default/nodes.dmp"
params.refseqkaijunames = "/storage/aku048/database/refseq/kaiju1.7.3/refseq_arch_bact_default/names.dmp"
params.refseqkraken = "/storage/aku048/database/refseq/kraken2/"
params.mar_refnucl = "/storage/aku048/database/mar5/kraken_mar_ref_nucl/database/"
params.mar_refprot = "/storage/aku048/database/mar5/kaiju_mar_ref_protein/marref_proteins_V5.fmi"
params.mar_dbnucl = "/storage/aku048/database/mar5/kraken_mar_db_nucl/database/"
params.mar_dbprot = "/storage/aku048/database/mar5/kaiju_mar_db_protein/masked_segmar5_db_prot.fmi"
params.mar_ref_dbnucl = "/storage/aku048/database/mar5/kraken_mar_ref_db_nucl/database/"
params.mar_ref_dbprot = "/storage/aku048/database/mar5/kaiju_mar_ref_db_protein/mar5_ref_db_prot.fmi"
*/

println """\
         M A R I N E - METAGENOMICS P I P E L I N E
         ==========================================
         reads_illumina   : ${params.illumina}
         reads_bgi        : ${params.bgi}
         outdir_illumina  : ${params.outputill}
         outdir_bgi       : ${params.outputbgi}
         """
         .stripIndent()



// Include modules for each tool
include { First_Fastqc; Optical; Clumpify; Fastqscreen; Repair; Second_Fastqc; Multiqc;  } from './processes/preprocessing_illumina.nf' 
include { First_Fastqcb; Clumpifyb; Repairb; Second_Fastqcb; Multiqcb; } from './processes/preprocessing_bgi.nf'
//include { Kaiju; Kaiju as Kaiju_bgi} from './processes/taxonomy_illumina.nf'
//include { Kaiju_refseq; Kaiju_refseq as Kaiju_refseq_bgi} from './processes/taxonomy_kaiju_refseq.nf'
//include { Kaiju_name} from './processes/names_taxonomy_kaiju.nf'

//include { Kraken2_marrefdb; Kraken2_marrefdb as Kraken2_marrefdb_bgi} from './processes/taxonomy_kraken_marref.nf'
//include { Kraken2_refseqdb; Kraken2_refseqdb as Kraken2_refseqdb_bgi} from './processes/taxonomy_kraken_refseq.nf'
//include {Assembly; Assembly as Assembly_bgi} from './processes/assembly.nf'
//include {Mapping; Mapping as Mapping_bgi; Maxbin; Metabat} from './processes/binning.nf'
//include {AddTaxonNames_kaiju} from './processes/phyloseq.nf'


//Preprocessing

//set the reads channel
Channel.fromFilePairs( "${params.illumina}/*_{fw,rv}.fastq.gz",  checkExists:true )
    .set{ illumina_read_pairs_ch }
Channel.fromFilePairs( "${params.bgi}/*_{1,2}.fastq.gz",  checkExists:true )
    .set{ bgi_read_pairs_ch }
    
//Run the main workflow    
workflow illumina_pre {
  main:
    First_Fastqc(illumina_read_pairs_ch ) 
  /*  Optical(illumina_read_pairs_ch)
    Clumpify(Optical.out.optical, params.adapter_ill)
    Fastqscreen(Clumpify.out.trim_r1, params.fastqscreen) | Repair | Second_Fastqc 
    Multiqc(First_Fastqc.out.collect(), Second_Fastqc.out.collect())
  emit:
    data = Repair.out */
}


workflow bgi_pre {
  main:
    First_Fastqcb(bgi_read_pairs_ch)
 /*   Clumpifyb(bgi_read_pairs_ch, params.adapter_bgi) 
    Repairb(Clumpifyb.out.trim) | Second_Fastqcb 
    Multiqcb(First_Fastqcb.out.collect(), Second_Fastqcb.out.collect())
  emit:
    data = Repairb.out  */
}

/*
//Taxonomic classification
workflow taxonomy_ill{
    take: my_data
    main:
      Kaiju(params.refseqkaijunodes, params.mar_ref_dbprot, my_data)
      Kaiju_refseq(params.refseqkaijunodes, params.refseqkaiju, my_data)
	  Kraken2_marrefdb(params.mar_ref_dbnucl, my_data)
	  Kraken2_refseqdb(params.refseqkraken, my_data)
    emit: 
      data = Kaiju.out
      data_Kaiju_refseq = Kaiju_refseq.out
	  data_Kraken2_marrefdb = Kraken2_marrefdb.out
	  data_Kraken2_refseqdb = Kraken2_refseqdb.out
}

workflow taxonomy_bgi{
    take: data
    main:
      Kaiju_bgi(params.refseqkaijunodes, params.mar_ref_dbprot, data)
      Kaiju_refseq_bgi(params.refseqkaijunodes, params.refseqkaiju, data)
	  Kraken2_marrefdb_bgi(params.mar_ref_dbnucl, data)
	  Kraken2_refseqdb_bgi(params.refseqkraken, data)
	emit: 
      data_Kaiju_bgi = Kaiju_bgi.out
      data_Kaiju_refseq_bgi = Kaiju_refseq_bgi.out
	  data_Kraken2_marrefdb_bgi = Kraken2_marrefdb_bgi.out
	  data_Kraken2_refseqdb_bgi = Kraken2_refseqdb_bgi.out
}


//Comperative Taxonomics analysis
workflow comperative_taxon_kaiju{
	take: ill_mar
	main:
		Kaiju_name(params.refseqkaijunodes, params.refseqkaijunames, ill_mar)
	emit:
		name_Kaiju = Kaiju_name.out
}



//Assembly
workflow assembly_ill{
    take: my_data
    main:
      Assembly(my_data)
    emit:
      contigs = Assembly.out
}

workflow assembly_bgi{
    take: data
    main:
      Assembly_bgi(data)
    emit:
      contig = Assembly_bgi.out
}


//Binning
workflow binning_ill{
    take: contigs
          data
    main:
      Mapping(contigs, data)
      Maxbin(contigs, Mapping.out.abun, data)
      Metabat(contigs, Mapping.out.sort, data)
}

workflow binning_bgi{
    take: my_data
          data1
    main:
      Mapping_bgi(my_data, data1)
}

*/

//Final Call
workflow {
    
    main:
      illumina_pre()
      bgi_pre()
      
/*      
      taxonomy_ill(illumina.out.data)  
      assembly_ill(illumina.out.data)
      binning_ill(assembly_ill.out.contigs, illumina.out.data)
      

      taxonomy_bgi(bgi.out.data)      
      assembly_bgi(bgi.out.data)
      binning_bgi(assembly_bgi.out.contig, bgi.out.data)

      comperative_taxon_kaiju(taxonomy_ill.out[0])	
		  comperative_taxon_kraken(taxonomy_ill.mar, taxonomy_ill.refseq, taxonomy_bgi.mar, taxonomy_bgi.refseq) , taxonomy_ill.out.Kaiju_refseq, taxonomy_bgi.data_Kaiju_bgi, taxonomy_bgi.data_Kaiju_refseq_bgi
*/
}


workflow.onComplete {
	log.info ( workflow.success ? "\nAnalysis Done!\n" : "Oops ... something went wrong" )}
 

//Done
