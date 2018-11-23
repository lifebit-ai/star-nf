params.star_index = "s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg19/Sequence/STARIndex/"

params.project = "SRP033351"

params.resultdir = 'results'

projectSRId = params.project

star_index = Channel.fromPath(params.star_index)
					.ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }

int threads = Runtime.getRuntime().availableProcessors()

process getSRAIDs {
	
	cpus 1

	input:
	val projectID from projectSRId
	
	output:
	file 'sra.txt' into sraIDs
	
	script:
	"""
	esearch -db sra -query $projectID  | efetch --format runinfo | grep SRR | cut -d ',' -f 1 > sra.txt
	"""
}

sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }

process fastqDump {

	publishDir params.resultdir, mode: 'copy'

	cpus threads

	input:
	val id from singleSRAId

	output:
	file '*.fastq.gz' into reads

	script:
	"""
	parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip
	"""	
}

process star {

	publishDir params.resultdir, mode: 'copy'

	cpus threads

	input:
	file read from reads
	file index from star_index.collect()

	output:
	file '*.bam' into alignedReads

	script:
	readName = read.toString() - ~/(\.fastq\.gz)?$/
	
	"""
	STAR --genomeDir $index --readFilesIn $read --runThreadN ${task.cpus} --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $readName

	"""
}

process count {

	publishDir params.resultdir, mode: 'copy'

	cpus threads

	input:
	file('*') from alignedReads.collect()

	output:
	file('counts.RData') into countData

	script:
	"""
        #!/usr/bin/env Rscript

	library(Rsubread)
	allbams = list.files(pattern=\"*.bam\")
	hg19ann <- getInBuiltAnnotation(\"hg19\")
	counts <- featureCounts(files=allbams, annot.ext=hg19ann, nthreads=${task.cpus})
	colnames(counts\$counts) <- allbams
	save(counts, file=\"counts.RData\")
	"""
}
