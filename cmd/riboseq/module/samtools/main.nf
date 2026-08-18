nextflow.enable.dsl=2

process index_samtools {
	tag "index_samtools"
	label "index_samtools"
	publishDir "$params.outdir/results/mapping_genome/", pattern:"*.bai", mode: 'copy', saveAs: { filename -> "${pair_ID}.bai"} 
	
	input:
		tuple val(pair_ID), path (file)
	
	output:
		tuple val(pair_ID), path ("*.bam*"), emit: index_bai_ch

		
	script:
	"""
		samtools index $file \
		$params.option_samtools
	"""
}
