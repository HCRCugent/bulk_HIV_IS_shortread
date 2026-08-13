# bulk_HIV_IS_shortread

The pipeline requires the following which are not included in this set:

	-	UCSC blat command line tool.
	-	hg38 genome as one fasta file per chromosome. The files 
		should be named as chr1.fa, chr2.fa...chr22.fa, chrX.fa, chrY.fa, chrM.fa, chrUn(*).fa.
  - usearch package (tested with: usearch11.0 (1).667_i86linux64)
	
		
The files included in this set are:
	
	-	1_demultiplex_inlinePE1.pl
  -	1_demultiplex_inlinePE1_PE2.pl
	-	2_blat_pipeline_usearch.pl 
	-	IlluminaBarcodes.pm
	-	ISAbatchFilter.pm
	-	ISAblatFilter.pm
	-	ISAdataReformat.pm
	- 11.ooc (required for Blat)
	-	example_settings_1_demultiplex_inlinePE1.txt
  -	example_settings_1_demultiplex_inlinePE1_PE2.txt
	-	README
