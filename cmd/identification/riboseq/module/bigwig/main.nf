nextflow.enable.dsl=2

process bigwig {
	tag "bigwig"
	label "bigwig"
	publishDir "$params.outdir/results/bigwig_genome/", pattern:"*.bw", mode: 'copy', saveAs: { filename -> "${pair_ID}.bw"} 
	
	input:

	tuple val(pair_ID), file (bam)
	file (idx)
		
	output:
		tuple val(pair_ID), file("*.bw"), emit : bw_ch
		
	script:
		
	"""
	
		bamCoverage \
		--binSize 1 \
		--bam ${bam} \
		--outFileName ${pair_ID}.bw \
		--verbose \
		--numberOfProcessors ${task.cpus} \
		$params.option_bigwig

	"""
}
