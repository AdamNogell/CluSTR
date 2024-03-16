# Main script:
Arguments:

	$1 = reads1.fq.gz
	$2 = reads2.fq.gz
	$3 = primers.txt
	$4 = adapters.txt
	$5 = reference_genome.fasta
primers.txt ($3) needs to be in format: 

  	FW 	Sequence
	RV  	Sequence
			 
adapters.txt ($4) needs to be in format:

	5'Adapter
	3'Adapter
