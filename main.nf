#!/usr/bin/env nextflow

nextflow.enable.dsl = 2 

// Make parameter list then move to yalm or config file

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:
      nextflow run /hpcnfs/data/GN2/fgualdrini/tools/NEXTFLOW_PIPES/nfSQuIRE --input design.csv --build mm10
    Mandatory arguments:
      --fastq_file_path [file]      Path to fastq.gz files
      --build [str]                 biuld version i.e. mm10 hg38 etc..
      --read_length         length of the reads
      --trim3               trim bases

    """.stripIndent()
}


// Show help message
if (params.help) {
    helpMessage()
    exit 0
}



workflow {

    if (params.fastq_file_path) { 
        ch_fastq_file_path = Channel.fromFilePairs( params.fastq_file_path + '*.R{1,2}.fastq.gz',checkIfExists:true) 
        } else { 
            exit 1, 'Specify the path to PE fastqs' 
            }

    if (params.build) { ch_build = params.build } else { exit 1, 'build must be one of the goldenPath UCSC genomes i.e. mm10, hg38 etc..' }
    if (params.read_length) { ch_read_length= params.read_length } else { exit 1, 'Specify the read length --read_length' }
    if (params.trim3) { ch_trim3 = params.trim3 } else { ch_trim3 = 0 }
    if (params.strandedness) { ch_strandedness = params.strandedness } else { ch_strandedness = 0 }
    if (params.em) { ch_em = params.em } else { ch_em = 'auto' }

    // Check input path parameters to see if they exist
    checkPathParamList = [params.fetch_folder, params.star_index,params.clean_folder]
    for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

    ch_fetch_folder = Channel.empty()
    ch_star_index = Channel.empty()

    if (!params.fetch_folder) {
        // SQuIRE FETCH files from golden path
        FETCH( ch_build )
        ch_fetch_folder = FETCH.out.squire_fetch
    } else {
        ch_fetch_folder = params.fetch_folder
    }

    if (!params.star_index) {
        // // Index with STAR
        STAR_INDEX( ch_build,
                    ch_fetch_folder)
        ch_star_index = STAR_INDEX.out.index
    } else {
        ch_star_index = params.star_index
    }


    ch_clean_folder = Channel.empty()
    if (!params.clean_folder ) {
        // // SQuIRE Clean
        CLEAN(  ch_build,
                ch_fetch_folder )
        ch_clean_folder = CLEAN.out.squire_clean
    } else {
        ch_clean_folder = params.clean_folder
    }

    // // SQuIRE Map
    MAP(    ch_fetch_folder,
            ch_star_index,
            ch_read_length,
            ch_trim3,
            ch_fastq_file_path
            )

    ch_mapped_bam = MAP.out.mapped_bam
    
    // Squire Count
    COUNT(
            ch_mapped_bam,
            ch_clean_folder,
            ch_fetch_folder,
            ch_read_length,
            ch_build,
            ch_strandedness,
            ch_em
    )

    // Squire Draw - generate BigWig and BedGraph
    DRAW(
            ch_mapped_bam,
            ch_fetch_folder,
            ch_build,
            ch_strandedness
    )

    // There are three sets of files to combine:
    // _subFcounts.txt
    // _TEcounts.txt
    // _refGenecounts.txt

    // Run with Rscript in bin - collect the individual files given a list of path - samples IDs are in .../squire_count_*/
    // append column wide per sample - In input only the path of the counts.

    
    COUNT.out.squire_count
            .collect{it[1]}
            .set { ch_count_all_path }

    NORMDEG (

        ch_count_all_path

    )


}

// Use the Fetch 

process FETCH {
    
    label 'process_high'
    container 'docker://lizatym/squire'
    echo true

    input:
    val build

    output:
    path squire_fetch       , emit: squire_fetch

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    squire Fetch \\
        -b $build \\
        -f \\
        -c \\
        -r \\
        -g \\
        -p $task.cpus \\
        -k
    """
}

process STAR_INDEX {

    label 'process_high'
    container 'docker://lizatym/squire'
    echo true


    input:
    val build
    path squire_fetch

    output:
    path "${build}_STAR"       , emit: index

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    def memory = task.memory ? "--limitGenomeGenerateRAM ${task.memory.toBytes() - 500000000}" : ''
    """
        mkdir ${build}_STAR
        STAR \\
            --runMode genomeGenerate \\
            --genomeDir ${build}_STAR/ \\
            --genomeFastaFiles ${squire_fetch}/${build}.chromFa/*.fa \\
            --runThreadN $task.cpus \\
            $memory \\
            $args
    """
}

process CLEAN {

    label 'process_medium'
    container 'docker://lizatym/squire'
    echo true



    input:
    val build
    path squire_fetch
    
    output:
    path squire_clean       , emit: squire_clean

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    mkdir squire_clean
    squire Clean \\
     -r $squire_fetch/${build}_rmsk.txt \\
     -b $build \\
     -o squire_clean/
    """
}


process MAP {

    label 'process_high'
    container 'docker://lizatym/squire'
    echo true


    input:
    path squire_fetch
    path index
    val read_length
    val trim3
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('mapped_*')       , emit: mapped_bam  // mapped is a folder!

    when:
    task.ext.when == null || task.ext.when
    
    script:

    def avail_mem = task.memory ? ((task.memory.toGiga() - 6) / task.cpus).trunc() : false
    def sort_mem = avail_mem && avail_mem > 2 ? '-m ${avail_mem}G' : '-m 1G'

    """
    mkdir mapped_${sample_id}
    STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${reads} \\
        --runThreadN ${task.cpus} \\
        --sjdbGTFfile ${squire_fetch}/*_refGene.gtf \\
        --sjdbOverhang \$((${read_length} -1)) \\
        --twopassMode Basic \\
        --readFilesCommand gunzip -c \\
        --outFilterMultimapNmax 100 \\
        --winAnchorMultimapNmax 100 \\
        --alignEndsType EndToEnd \\
        --alignEndsProtrude 100 DiscordantPair \\
        --outFilterScoreMinOverLread 0.4 \\
        --outFilterMatchNminOverLread 0.4 \\
        --chimSegmentMin ${read_length} \\
        --outFileNamePrefix ./mapped_${sample_id}/${sample_id}. \\
        --outSAMtype BAM Unsorted \\
        --outSAMattributes All \\
        --outSAMstrandField intronMotif \\
        --outSAMattrIHstart 0 \\
        --clip3pNbases $trim3

     
    samtools sort -@ ${task.cpus} ${sort_mem} -o ./mapped_${sample_id}/${sample_id}.sort.out.bam -T ./mapped_${sample_id}/${sample_id}.sort.out ./mapped_${sample_id}/${sample_id}*.out.bam
    samtools index ./mapped_${sample_id}/${sample_id}.sort.out.bam
    rm -rf ./mapped_${sample_id}/*._STAR*
    rm -rf ./mapped_${sample_id}/*Chimeric*
    rm -rf ./mapped_${sample_id}/*Aligned.out.bam
    rm -rf ./mapped_${sample_id}/*.out.tab


    """
}

process COUNT {
    
    label 'process_medium'
    container 'docker://lizatym/squire'
    echo true

    input:
    tuple val(sample_id), path(mapped_bam) 
    path squire_clean
    path fetch_folder
    val read_length
    val build
    val strandedness
    val em

    output:
    tuple val(sample_id), path('squire_count_*')     , emit: squire_count

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    mkdir squire_count_${sample_id}
    squire Count \\
        -m ./${mapped_bam} \\
        -c ${squire_clean} \\
        -o ./squire_count_${sample_id} \\
        -t ./squire_count_${sample_id} \\
        -f ${fetch_folder} \\
        -r ${read_length} \\
        -b ${build} \\
        -p ${task.cpus} \\
        -s ${strandedness} \\
        -e ${em}  
    """
    
}


process DRAW {
    
    label 'process_high'
    container 'docker://lizatym/squire'
    echo true

    input:
    tuple val(sample_id), path(mapped_bam) 
    path fetch_folder
    val build
    val strandedness


    output:
    tuple val(sample_id), path('squire_draw_*')     , emit: squire_draw

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    mkdir squire_draw_${sample_id}
    squire Draw \\
            -f ${fetch_folder} \\
            -m ./${mapped_bam} \\
            -o ./squire_draw_${sample_id} \\
            -n ${sample_id} \\
            -s ${strandedness} \\
            -b ${build} \\
            -p ${task.cpus} 
    """
    
}


process NORMDEG {
    
    label 'process_high'
    container 'docker://lizatym/squire'
    echo true

    input:
    path ch_count_all_path
    val build
    val strandedness


    output:
    tuple val(sample_id), path('squire_draw_*')     , emit: squire_draw

    when:
    task.ext.when == null || task.ext.when
    
    script:
    def args = task.ext.args ?: ''
    """
    mkdir squire_draw_${sample_id}
    squire Draw \\
            -f ${fetch_folder} \\
            -m ./${mapped_bam} \\
            -o ./squire_draw_${sample_id} \\
            -n ${sample_id} \\
            -s ${strandedness} \\
            -b ${build} \\
            -p ${task.cpus} 
    """
    
}
