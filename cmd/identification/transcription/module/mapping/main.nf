nextflow.enable.dsl=2


/*
 * Process 1 : construction du génome indexé par l'outil STAR
 */  

process indexation_genome_STAR {
	label "indexation_genome_STAR"
    tag "${indexation_genome_STAR}"
    
  input:
      path(genome)
       
  output:
      path results_index ,emit : index


  script:
  """
    mkdir results_index 
    
    STAR \
	--runThreadN ${task.cpus} \
	--genomeDir results_index \
	--runMode genomeGenerate \
	--genomeFastaFiles $genome \
	$params.option_star_index
	
  """
   }


/*
 * Process 2 : mapping des reads sur le génome de référence indexé
 */
 
process mapping_STAR {
	label "mapping_STAR"
	tag "${mapping_STAR}"
	publishDir "$params.outdir/results/mapping_genome/", pattern:"*sortedByCoord.out.bam", mode: 'copy', saveAs: { filename -> "${pair_ID}.bam" }
	publishDir "$params.outdir/results/mapping_log/", pattern:"*Log.final.out", mode: 'copy', saveAs: { filename -> "${pair_ID}.out" }
	publishDir "$params.outdir/results/mapping_transcriptome/", pattern:"*toTranscriptome.out.bam", mode: 'copy', saveAs: { filename -> "${pair_ID}.transcriptome.bam" }
	
	
	input:
		path(results_index)
		tuple val(pair_ID), path(reads)
		path(annotation)
	
	output:
		file "*"
		tuple val(pair_ID), path ("${pair_ID}Aligned.sortedByCoord.out.bam") , emit: bam
		tuple val (pair_ID), path ("${pair_ID}Aligned.toTranscriptome.out.bam"), emit: transcri_bam

	script:
	"""	
		STAR \
		--runThreadN ${task.cpus} \
		--readFilesIn $reads \
		--genomeDir $results_index \
		--outFileNamePrefix $pair_ID \
		--outSAMtype BAM SortedByCoordinate \
		--quantMode TranscriptomeSAM \
		--sjdbGTFtagExonParentTranscript Parent \
		--sjdbGTFfile $annotation \
		--alignEndsType EndToEnd \
		--twopassMode  Basic \
		$params.option_star_mapping
	 """
}
