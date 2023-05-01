/* 
 * pipeline input parameters 
 */

params.samplesheet = "${baseDir}/Samplesheets/samples_test.csv"
params.outdir = "${baseDir}/results"
params.index = "/Zulu/bnolan/Indexes/STARindexhg38/"
params.gtf = "/Zulu/bnolan/genome/human/gtf/hg38.ensGene.gtf"
params.threads = "4"
params.notrim = ""
params.singleend = ""


log.info """\
         ===================================
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         outdir       : ${params.outdir}
         samplesheet  : ${params.samplesheet}
         threads      : ${params.threads}
         index        : ${params.index}
         gtf          : ${params.gtf}
         single-end   : ${params.singleend}

         trim         : ${params.notrim}

         """
         .stripIndent()


// Parse samplesheet and create reads channel; 
// dependent on single-end or paired-end sequencing   

if(params.singleend){
    Channel
        .from ( file(params.samplesheet) )
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true) ]] }
        .set { read_pairs_ch }
} else {
    Channel
        .from ( file(params.samplesheet) )
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.sample_id, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]] }
        .set { read_pairs_ch }
}



// Create channel for Index
Channel
        .value(file(params.index))
        .set{index_ch}

// Create channel for gtf
Channel
        .value(file(params.gtf))
        .set{gtf_ch}


// // Concatenate reads of same sample [e.g. additional lanes, resequencing of old samples]
//nfcore
read_pairs_ch
        .groupTuple()
        .branch {
            meta, fastq ->
                single  : fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastq }  


// Concatenate 
process CAT_FASTQ {
    //nf-core
    tag "$sample_id"
    publishDir "${params.outdir}/fastqMerged", mode: 'copy'

    input:
    tuple val(sample_id), path(reads, stageAs: "input*/*") from ch_fastq.multiple

    output:
    tuple val(sample_id), path("*.merged.fastq.gz") into cat_out_ch


    script:
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]

    if (params.singleend) {
            if (readList.size >= 1) {
                """
                cat ${readList.join(' ')} > ${sample_id}.merged.fastq.gz
                """
            }
        } else {
            if (readList.size >= 2) {
                def read1 = []
                def read2 = []
                readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
                """
                cat ${read1.join(' ')} > ${sample_id}_1.merged.fastq.gz
                cat ${read2.join(' ')} > ${sample_id}_2.merged.fastq.gz
                """
            }
        }
}

//mix merged reads with singles
cat_out_ch
        .mix(ch_fastq.single)
        .into { cat_merged_ch; cat_merged_ch2; cat_merged_ch3 }



//  Run fastQC to check quality of reads files
process fastqc {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc", pattern:"{*.html,fastqc_${sample_id}_logs}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads) from cat_merged_ch

    output:
    path("fastqc_${sample_id}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}

// Trimming reads with Trim Galore

process trimming {
    tag "$key"

    publishDir "${params.outdir}/trimmed", pattern: "*.fq.gz", mode: 'copy'

    input: 
    tuple val(key), path(reads) from cat_merged_ch2

    when:
    !params.notrim

    output:
    tuple val(key), path("*.fq.gz") into reads_trimmed_ch

    script:

    if(params.singleend){

        """
        trim_galore \\
                    $reads \\
                    --basename $key \\
                    --cores 1 
        """

    } else {

        """
        trim_galore \\
                    --paired \\
                    ${reads[0]} \\
                    ${reads[1]} \\
                    --basename $key \\
                    --cores 1 
        """
    }

    
}

reads_out_ch = params.notrim ? cat_merged_ch3 : reads_trimmed_ch


// Align reads to index 

process STAR {
    tag "$key"
    publishDir "${params.outdir}/STAR/$key/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(reads) from reads_out_ch
    path(index) from index_ch  

    output:
    tuple val(key), path("*bam") into bam_ch, bam_ch2
    path("*Log.final.out") into align_report_ch 

    script:

        """
        STAR \\
                --runThreadN $params.threads \\
                --genomeDir $index \\
                --readFilesIn $reads \\
                --readFilesCommand zcat \\
                --outFileNamePrefix $key \\
                --outStd SAM \\
                | samtools view -@ $params.threads -bhS -o ${key}.bam -
                
        """  
}

 

// Samtools stats for summary statistics on bwa-meth alignment
process samtools_stat_flagstat {
    tag "$key"
    publishDir "${params.outdir}/samtools/stats/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_ch


    output:
    path('*stats') into samtools_stats_ch
    path('*flagstat') into samtools_flag_ch


    script:
    """
    samtools stats --threads ${params.threads} \\
                   $bam \\
                   > ${key}.stats

    samtools flagstat \\
                        --threads ${params.threads} \\
                        $bam \\
                        > ${key}.flagstat

    """
}



process samtools_sort {
    tag "$key"
    publishDir "${params.outdir}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_ch2
    
    output:
    tuple val(key), path('*sort.bam') into bam_sorted_ch, bam_sorted_ch2, bam_sorted_ch3

    script:
    """
    samtools sort $bam > ${key}.sort.bam
    """
} 

process samtools_index {
    tag "$key"
    publishDir "${params.outdir}/samtools/bam/${key}/", pattern:'*', mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sorted_ch

    output:
    tuple val(key), path('*bai') into bam_indexed_ch

    script:

    """
    samtools index $bam
    """
}


// Add bai to bam channel for bamCoverage
bam_sorted_ch2
            .join(bam_indexed_ch)
            .set{bam_sorted_indexed_ch}



process bamCoverage {
    tag "$key"
    publishDir "${params.outdir}/bamCoverage/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam), path(bai) from bam_sorted_indexed_ch

    output:
    tuple val(key), path("*.bw") into bw_ch

    script:
    """
    bamCoverage -b $bam -o ${key}.bw --normalizeUsing BPM
    """
}


// // Count transcripts 
process stringtie {
    tag "$key"
    publishDir "${params.outdir}/stringtie/$key", pattern:"*", mode: 'copy'

    input:
    tuple val(key), path(bam) from bam_sorted_ch3
    path(gtf) from gtf_ch


    output:
    tuple val(key), path("*gtf") into stringtie_ch
    tuple val(key), path("*.ctab") into ctab_ch
    tuple val(key), path("*.tab") into tab_ch

    script:

    """
    stringtie \\
                $bam \\
                -o ${key}.gtf \\
                -p $params.threads \\
                -G $gtf \\
                -l $key \\
                -eB \\
                -A ${key}_gene_abund.tab 
    """
}




// Create multiqc report channel
  multiqc_ch = fastqc_ch
                    .mix(align_report_ch)                    
                    .mix(samtools_flag_ch)
                    .mix(samtools_stats_ch)
                    .collect()


process multiqc {

    publishDir "${params.outdir}/multiqc/", pattern:"multiqc_report.html", mode:'copy'
       
    input:
    path('*') from multiqc_ch
    
    output:
    path('multiqc_report.html')  
     
    script:
    """
    multiqc . 
    """
} 


workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc/multiqc_report.html\n" : "Oops .. something went wrong" )
}