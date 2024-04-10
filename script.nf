#!/usr/bin/env nextflow

//help function for the tool
def show_help (){
  log.info IARC_Header()
  log.info tool_header()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run iarcbioinfo/TCR-BCR-nf -singularity [OPTIONS]

    Mandatory arguments:
      --input_folder         [file] Folder with bam files to be processed
      --IMGTC_fasta          [file] fasta file of reference genome from the international ImMunoGeneTics information system [file human_IMGT+C.fa from https://github.com/liulab-dfci/TRUST4]
      --bcrtcr_fasta         [file] fasta file with reference BCR and TCR regions [see https://github.com/liulab-dfci/TRUST4 for generation]
    Optional arguments:
      --output_folder        [string] name of output folder
      --cpu                  [Integer]  Number of CPUs[def:2]
      --mem 		         [Integer] Max memory [def:8Gb]
      --barcode              [flag] Run trust4 using a specific barcode (for single cell data only, usually BC for CellRanger output)
      """.stripIndent()
}


//run trust4
process trust {
  cpus params.cpu
  memory params.mem+'G'
  tag { bam }

  publishDir params.output_folder, mode: 'copy'
  input:
  path(bam)
  path(bcrtcr_fasta)
  path(IMGTC_fasta)
  val(singlecell)

  output:
  tuple val(tumor_id), file("report-${tumor_id}-hla.json")
  script:
       """
       run-trust4 -b {!bam}  -f !{bcrtcr_fasta} --ref !{IMGTC_fasta} !{singlecell}
       """
}


// DSL2 workflow to run the processes
workflow{
  //display help information
  if (params.help){ show_help(); exit 0;}
  //display the header of the tool
  log.info IARC_Header()
  log.info tool_header()
  //Check mandatory parameters
  assert (params.input_folder != null) : "please specify --bam file"
  assert (params.bcrtcr_fasta != null ) : "please specify --bcrtcr_fasta"
  assert (params.IMGTC_fasta != null ) : "please specify --IMGTC_fasta"


  //channels for reference genome
  IMGTC_fasta = Channel.value(file(params.IMGTC_fasta)).ifEmpty{exit 1, "reference file not found: ${params.IMGTC_fasta}"}
  bcrtcr_fasta = Channel.value(file(params.bcrtcr_fasta)).ifEmpty{exit 1, "reference file not found: ${params.bcrtcr_fasta}"}
  
  //channel with bam files
  bams = Channel.fromPath(params.input_folder+"/*.bam")

  //to add barcode ID; usually BC
  if(params.barcode){
    singlecell="--barcode "+params.barcode
  }else{
    singlecell=" "
  }

  print_params()

  //run TCR and BCR genotyping
  trust(bams,bcrtcr_fasta,IMGTC_fasta,singlecell)
}


// print the calling parameter to the log and a log file
def print_params () {
  //software versions
  def software_versions = ['trust4'   : '1.1.0']
  //we print the parameters
  log.info "\n"
  log.info "-\033[2m------------------Calling PARAMETERS--------------------\033[0m-"
  log.info params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"
  log.info "-\033[2m------------------Software versions--------------------\033[0m-"
  log.info software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
  log.info "-\033[2m--------------------------------------------------------\033[0m-"
  log.info "\n"


  //we print the parameters to a log file
   def output_d = new File("${params.output_folder}/nf-pipeline_info/")
   if (!output_d.exists()) {
       output_d.mkdirs()
   }
   def output_tf = new File(output_d, "run_parameters_report.txt")
   def  report_params="------------------Calling PARAMETERS--------------------\n"
        report_params+= params.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="\n------------------NEXTFLOW Metadata--------------------\n"
        report_params+="nextflow version : "+nextflow.version+"\n"
        report_params+="nextflow build   : "+nextflow.build+"\n"
        report_params+="Command line     : \n"+workflow.commandLine.split(" ").join(" \\\n")
        report_params+="\n--------------------------------------------------------\n"
        report_params+="-----------------Software versions--------------------\n"
        report_params+=software_versions.collect{ k,v -> "${k.padRight(18)}: $v"}.join("\n")
        report_params+="\n--------------------------------------------------------\n"

   output_tf.withWriter { w -> w << report_params}
}


//this use ANSI colors to make a short tool description
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        TCR-BCR-nf: Pipeline to genotype Tcell and Bcell receptors from RNA-seq data (${workflow.manifest.version})
        """
}

//header for the IARC tools
// the logo was generated using the following page
// http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}