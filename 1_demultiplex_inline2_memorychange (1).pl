#!/usr/bin/perl

# 1_demultiplex_inline2.pl
# Dr. Xiaolin Wu & Daria W. Wells - update Dr. Jolien Vermeire

# This script is used to filter Illumina paired-end ISA sequence data from the MiSeq, 
# HiSeq, or NextSeq. It will search for the user-specified barcodes, remove unspecified barcodes, and 
# and split the reads at the LTR/Linker junctions for downstream processing. It 
# requires knowledge of the bases immediately upstream to the junction. 
#
# Optional PE2 inline index support:
#   If <PE2barcode_present> is set to TRUE in the settings file, the script will also
#   extract and match an 8nt inline index from the start of Read 2 against the lists
#   of PE2 barcodes. When enabled:
#     - The UMI is extracted from bases 8-17 of Read 2 (instead of 0-9).
#     - <PE1barcodes> and <PE2barcodes> hold THIS sample's expected barcodes
#       (parallel lists of equal length; usually length 1).
#     - <PE1barcodes_other> and <PE2barcodes_other> (optional) hold barcodes used by
#       OTHER samples on the same flowcell. Reads carrying these are likely the result
#       of index hopping. The two "other" lists do not need to be the same length.
#     - Reads are sorted into:
#         * Both expected                  -> <LTR>/PE1_<exp>_PE2_<exp>.txt
#         * Any combo containing an "other" barcode in either or both reads
#                                          -> <LTR>/mismatched_reads/PE1_<X>_PE2_<Y>.txt
#                                             (one file per specific (X,Y) combination)
#         * Only PE1 matched (expected)    -> <LTR>/single_index/PE1_expected.txt
#         * Only PE1 matched (other)       -> <LTR>/single_index/PE1_other.txt
#         * Only PE2 matched (expected)    -> <LTR>/single_index/PE2_expected.txt
#         * Only PE2 matched (other)       -> <LTR>/single_index/PE2_other.txt
#     - Mismatched pairs are reported per-combination. Single-index reads are reported
#       split by expected vs. other.
#   If <PE2barcode_present> is FALSE or absent, the script behaves exactly as before:
#     - The UMI is extracted from bases 0-9 of Read 2.
#     - <PE2barcodes> is used only for output file naming (first entry only).
#     - The "_other" barcode lists, if present, are ignored.

use strict;

#use IO::Zlib qw(:gzip_external 1);
use Carp;
$| = 1;   # autoflush STDOUT and STDERR

# Helper: safe print that dies on write failure (catches out-of-disk, OOM partial writes, etc.)
sub safe_print {
	my ($fh, $data, $path) = @_;
	print $fh $data or die "Write to $path failed: $!\n";
}

use lib '/data/gent/vo/001/gvo00120/bulkISS/bulk_IS_illumina_v2_inlineindexPE2/'; # The path of the Illumina barcodes directory.
use IlluminaBarcodes; # The module that returns a hash of Illumina barcodes. Update as needed.

# Require settings file on the command line.
croak "No settings file specified on the command line." unless @ARGV;


# Import your settings from the settings file. Each value is a list of the settings under
# header, allowing you to look for multiple barcodes and analyze both LTRs or only one.
my (%settings, $header); # Initialize the settings hash and the header variable.
open my $settings_fh, "<$ARGV[0]" or die "Cannot read the settings file $ARGV[0]: $!"; # Open the settings file for reading.
# add whitespace trimming for settings files created on Windows
while (my $settings_line = <$settings_fh>) { # 
	chomp($settings_line);
	$settings_line =~ s/\r$//;                  # strip Windows CR if present
	if ($settings_line =~ /<([^>]*)>/) {        # Text surrounded by < and > indicates a settings type header, e.g. <LTR>.
		$header = $1; # Assign the matched text to $header.
		$header =~ s/^\s+|\s+$//g;              # trim header whitespace
		next;
	}
	next if $settings_line =~ /^\s*$/;          # skip blank lines
	$settings_line =~ s/^\s+|\s+$//g;            # trim leading/trailing whitespace
	push @{$settings{$header}}, $settings_line;
}
close $settings_fh; # Finished with settings file.

# Create hashes for to hold the barcodes and file names.
my %all_barcodes = IlluminaBarcodes::get_barcodes(); # Obtain the hash of all Illumina barcodes.

# Determine whether PE2 inline index matching is enabled.
my $use_PE2_inline;
if (defined $settings{PE2barcode_present} and
    uc($settings{PE2barcode_present}[0]) eq 'TRUE') {
	$use_PE2_inline = 1;
} else {
	$use_PE2_inline = 0;
}

# UMI offset depends on whether PE2 inline index is present at the start of Read 2.
my $UMI_offset;
if ($use_PE2_inline) {
	$UMI_offset = 8;
} else {
	$UMI_offset = 0;
}
my $UMI_length = 10;

# When PE2 inline matching is on, PE1barcodes and PE2barcodes must be parallel lists
# (same length); position N in each defines the expected pairing for sample N.
if ($use_PE2_inline) {
	my $n1 = scalar @{$settings{PE1barcodes}};
	my $n2 = scalar @{$settings{PE2barcodes}};
	croak "When PE2barcode_present=TRUE, <PE1barcodes> and <PE2barcodes> must have the same number of entries (got $n1 PE1 and $n2 PE2).\n"
		unless $n1 == $n2;
	croak "Sample list is empty.\n" if $n1 == 0;
}

# Hashes used for sample assignment / output file lookup:
#   %pe1_seq2name      : PE1 barcode sequence -> PE1 barcode name (covers expected and other)
#   %pe2_seq2name      : PE2 barcode sequence -> PE2 barcode name (covers expected and other)
#   %pe1_is_expected   : PE1 barcode sequence -> 1 if "expected" (this sample)
#   %pe2_is_expected   : PE2 barcode sequence -> 1 if "expected" (this sample)
#   %expected_pair     : PE1_seq -> expected partner PE2_seq (defines correct pairings)
#   %file_correct      : {LTR}{PE1_seq}                 -> filename for correct-pair reads
#   %file_mismatched   : {LTR}{PE1_seq}{PE2_seq}        -> filename for mismatched-pair reads
#   $file_pe1only_exp  : {LTR}                          -> filename for PE1-only (expected) reads
#   $file_pe1only_oth  : {LTR}                          -> filename for PE1-only (other)    reads
#   $file_pe2only_exp  : {LTR}                          -> filename for PE2-only (expected) reads
#   $file_pe2only_oth  : {LTR}                          -> filename for PE2-only (other)    reads
my (%pe1_seq2name, %pe2_seq2name,
    %pe1_is_expected, %pe2_is_expected,
    %expected_pair,
    %file_correct, %file_mismatched,
    %file_pe1only_exp, %file_pe1only_oth,
    %file_pe2only_exp, %file_pe2only_oth);

# Validate and resolve EXPECTED PE1 barcode names -> sequences.
foreach my $PE1_barcode (@{$settings{PE1barcodes}}) {
	croak "PE1 barcode name '$PE1_barcode' not found in IlluminaBarcodes.\n"
		unless exists $all_barcodes{$PE1_barcode};
	my $seq = $all_barcodes{$PE1_barcode};
	$pe1_seq2name{$seq} = $PE1_barcode;
	$pe1_is_expected{$seq} = 1;
}

# Validate and resolve EXPECTED PE2 barcode names -> sequences.
foreach my $PE2_barcode (@{$settings{PE2barcodes}}) {
	croak "PE2 barcode name '$PE2_barcode' not found in IlluminaBarcodes.\n"
		unless exists $all_barcodes{$PE2_barcode};
	my $seq = $all_barcodes{$PE2_barcode};
	$pe2_seq2name{$seq} = $PE2_barcode;
	$pe2_is_expected{$seq} = 1;
}

# Validate and resolve OTHER PE1 barcodes (optional).
if (defined $settings{PE1barcodes_other}) {
	foreach my $PE1_barcode (@{$settings{PE1barcodes_other}}) {
		croak "PE1 barcode name '$PE1_barcode' not found in IlluminaBarcodes.\n"
			unless exists $all_barcodes{$PE1_barcode};
		my $seq = $all_barcodes{$PE1_barcode};
		next if exists $pe1_is_expected{$seq};
		$pe1_seq2name{$seq} = $PE1_barcode;
	}
}

# Validate and resolve OTHER PE2 barcodes (optional).
if (defined $settings{PE2barcodes_other}) {
	foreach my $PE2_barcode (@{$settings{PE2barcodes_other}}) {
		croak "PE2 barcode name '$PE2_barcode' not found in IlluminaBarcodes.\n"
			unless exists $all_barcodes{$PE2_barcode};
		my $seq = $all_barcodes{$PE2_barcode};
		next if exists $pe2_is_expected{$seq};
		$pe2_seq2name{$seq} = $PE2_barcode;
	}
}

# Build expected-pair lookup using parallel-list semantics.
if ($use_PE2_inline) {
	for (my $i = 0; $i < @{$settings{PE1barcodes}}; $i++) {
		my $pe1_seq = $all_barcodes{$settings{PE1barcodes}[$i]};
		my $pe2_seq = $all_barcodes{$settings{PE2barcodes}[$i]};
		$expected_pair{$pe1_seq} = $pe2_seq;
	}
}

# Create file names to be used for data output
foreach my $LTR (@{$settings{'LTR'}}) {
	
	# Use the $LTR variable to create a file containing the primer and junction 
	# sequences for later use.
	mkdir $LTR unless -d $LTR;
	my $primer_file = "$LTR/${LTR}_primer.txt";
	my $primer_key  = "${LTR}_primer";
	my $junction_key = "${LTR}_junction";
	open my $primer_fh, '>', $primer_file
		or croak "Cannot open $primer_file for writing: $!";
	print $primer_fh @{$settings{$primer_key}}
		or croak "Cannot write to $primer_file: $!";
	print $primer_fh "\n@{$settings{$junction_key}}"
		or croak "Cannot write to $primer_file: $!";
	close $primer_fh
		or croak "Cannot close $primer_file: $!";
	
	if ($use_PE2_inline) {
		mkdir "$LTR/mismatched_reads";
		mkdir "$LTR/single_index";
		
		# Correct-pair files: one per sample (parallel-list pairing), placed in <LTR>/.
		for (my $i = 0; $i < @{$settings{PE1barcodes}}; $i++) {
			my $PE1_barcode = $settings{PE1barcodes}[$i];
			my $PE2_barcode = $settings{PE2barcodes}[$i];
			my $pe1_seq = $all_barcodes{$PE1_barcode};
			croak "Duplicate PE1 barcode '$PE1_barcode' (sequence $pe1_seq) in sample list for $LTR; ".
			      "each PE1 barcode may appear at most once per LTR.\n"
				if exists $file_correct{$LTR}{$pe1_seq};
			my $filename = "PE1_${PE1_barcode}_PE2_${PE2_barcode}.txt";
			croak "\nError: $LTR/$filename already exists\n" if -e "$LTR/$filename";
			$file_correct{$LTR}{$pe1_seq} = $filename;
		}
		
		# Mismatched-pair files: every (PE1, PE2) combination where at least one barcode
		# is "other", placed in <LTR>/mismatched_reads/.
		my @all_pe1_names = (@{$settings{PE1barcodes}},
		                     defined $settings{PE1barcodes_other} ? @{$settings{PE1barcodes_other}} : ());
		my @all_pe2_names = (@{$settings{PE2barcodes}},
		                     defined $settings{PE2barcodes_other} ? @{$settings{PE2barcodes_other}} : ());
		# Deduplicate name lists.
		my (%seen1, %seen2);
		@all_pe1_names = grep { !$seen1{$_}++ } @all_pe1_names;
		@all_pe2_names = grep { !$seen2{$_}++ } @all_pe2_names;

		foreach my $PE1_barcode (@all_pe1_names) {
			my $pe1_seq = $all_barcodes{$PE1_barcode};
			foreach my $PE2_barcode (@all_pe2_names) {
				my $pe2_seq = $all_barcodes{$PE2_barcode};
				my $pe1_exp = exists $pe1_is_expected{$pe1_seq};
				my $pe2_exp = exists $pe2_is_expected{$pe2_seq};
				next if $pe1_exp and $pe2_exp; # both expected -> not mismatched here
				my $filename = "PE1_${PE1_barcode}_PE2_${PE2_barcode}.txt";
				croak "\nError: $LTR/mismatched_reads/$filename already exists\n"
					if -e "$LTR/mismatched_reads/$filename";
				$file_mismatched{$LTR}{$pe1_seq}{$pe2_seq} = $filename;
			}
		}

		# Single-index files (single_index/): one combined file per category.
		foreach my $tag (qw(expected other)) {
			my $pe1_fname = "PE1_${tag}.txt";
			croak "\nError: $LTR/single_index/$pe1_fname already exists\n"
				if -e "$LTR/single_index/$pe1_fname";
			my $pe2_fname = "PE2_${tag}.txt";
			croak "\nError: $LTR/single_index/$pe2_fname already exists\n"
				if -e "$LTR/single_index/$pe2_fname";
		}
		$file_pe1only_exp{$LTR} = "PE1_expected.txt";
		$file_pe1only_oth{$LTR} = "PE1_other.txt";
		$file_pe2only_exp{$LTR} = "PE2_expected.txt";
		$file_pe2only_oth{$LTR} = "PE2_other.txt";

	} else {
		# Legacy: PE1_<name>_PE2_<single name>.txt directly under <LTR>/
		foreach my $PE1_barcode (@{$settings{PE1barcodes}}) {
			my $pe1_seq = $all_barcodes{$PE1_barcode};
			croak "Duplicate PE1 barcode '$PE1_barcode' (sequence $pe1_seq) in sample list for $LTR; ".
			      "each PE1 barcode may appear at most once per LTR.\n"
				if exists $file_correct{$LTR}{$pe1_seq};
			my $filename = "PE1_${PE1_barcode}_PE2_$settings{PE2barcodes}[0].txt";
			croak "\nError: $LTR/$filename already exists\n" if -e "$LTR/$filename";
			$file_correct{$LTR}{$pe1_seq} = $filename;
		}
	}
}


# Use subroutine choose_LTR_sub to return a subroutine reference to analyze either the 
# 3LTR only, 5LTR only, or both. Pass to choose_LTR_sub each LTR primer and junction 
# sequence, as well as the list of LTR's to be analyzed.
my $match_LTR_seq = choose_LTR_sub(
    ($settings{"3LTR_primer"}[0]   // ''),
    ($settings{"3LTR_junction"}[0] // ''),
    ($settings{"5LTR_primer"}[0]   // ''),
    ($settings{"5LTR_junction"}[0] // ''),
    @{$settings{"LTR"}},
);

# Create a list of the numbers from 1 to n in the format 001, ... 00n with 
# n being the number of Read 1 & Read 2 pairs for analysis.
my @files2read; # Initialize the array
for (my $i = 1; $i <= $settings{pairs}[0]; $i++) { # Perform this loop for all integers from 1 to n.
	my $pairs = sprintf("%03d", $i); # Format the number with leading zeros.
	push @files2read, $pairs; # Add the formatted number to @files2read.
}

# Create the file name template from the file name in the settings.
# The template is the text before R1 or R2, e.g. "DNAsample_S1_L001_"
my ($filenametemplate) =
    ($settings{file_name_template}[0] =~ /(.*)R\d(?:_\d{3})?\.fastq(?:\.gz)?$/)
    or croak "Incorrect FASTQ file name format. Check your settings file.\n";
	
# Open a file handle for printing the log file.
my $log_filename = "demultiplex_log.txt";
open my $log_fh, '>>', $log_filename
	or die "Unable to open $log_filename for appending: $!\n";
$log_fh->autoflush(1);   # flush every write so log stays current and ordered

# Print the file name template and the start time to the log file and/or terminal.
my $start_time = localtime;
print "\n$filenametemplate\n";
print $log_fh "Demultiplex script used: $0\n"; 
print $log_fh "Analysis started at $start_time\n\n";
print $log_fh "$filenametemplate\n";

# Add './' to the beginning of the file name template.
$filenametemplate="./$filenametemplate";

# Initialize the counter for the total number of paired end reads across all sequence files analyzed.
my $grand_total_reads = my $all_3LTR = my $all_5LTR = my $all_unidentified = 0;
my $all_correct = my $all_mismatched = 0;
my $all_pe1only_exp = my $all_pe1only_oth = my $all_pe2only_exp = my $all_pe2only_oth = 0;

# Per-combination mismatched-pair totals across all passes:
#   $mismatched_totals{$LTR}{$pe1_seq}{$pe2_seq} = N
my %mismatched_totals;

# Send the required variables and references to the extract_and_print_data subroutine.
foreach my $pairs (@files2read) {
	(my $total_reads, my $_3LTR_count, my $_5LTR_count, my $unidentified,
	 my $correct_count, my $mismatched_count,
	 my $pe1only_exp_count, my $pe1only_oth_count,
	 my $pe2only_exp_count, my $pe2only_oth_count,
	 my $mm_breakdown) = extract_and_print_data(
		$pairs, $filenametemplate,
		\%file_correct, \%file_mismatched,
		\%file_pe1only_exp, \%file_pe1only_oth,
		\%file_pe2only_exp, \%file_pe2only_oth,
		\%pe1_seq2name, \%pe2_seq2name,
		\%pe1_is_expected, \%pe2_is_expected,
		\%expected_pair,
		$settings{'LTR'}, $match_LTR_seq, $log_fh,
	    $use_PE2_inline, $UMI_offset, $UMI_length);
	$grand_total_reads += $total_reads;
	$all_3LTR           += $_3LTR_count;
	$all_5LTR           += $_5LTR_count;
	$all_unidentified   += $unidentified;
	$all_correct        += $correct_count;
	$all_mismatched     += $mismatched_count;
	$all_pe1only_exp    += $pe1only_exp_count;
	$all_pe1only_oth    += $pe1only_oth_count;
	$all_pe2only_exp    += $pe2only_exp_count;
	$all_pe2only_oth    += $pe2only_oth_count;
	# Aggregate per-combination mismatched totals across passes.
	foreach my $LTR (keys %$mm_breakdown) {
		foreach my $p1 (keys %{$mm_breakdown->{$LTR}}) {
			foreach my $p2 (keys %{$mm_breakdown->{$LTR}{$p1}}) {
				$mismatched_totals{$LTR}{$p1}{$p2} += $mm_breakdown->{$LTR}{$p1}{$p2};
			}
		}
	}
}

my %LTRs; # This hash is only used to determine which LTRs should be reported in the log
$LTRs{$_} = undef foreach @{$settings{'LTR'}};
# Print the number of paired end reads to the terminal and log file.

print "*" x 40 . "\n\n";
printf "Total 3LTR read pairs:\t%s\n", exists $LTRs{'3LTR'} ? $all_3LTR : "N/A";
printf "Total 5LTR read pairs:\t%s\n", exists $LTRs{'5LTR'} ? $all_5LTR : "N/A";
if ($use_PE2_inline) {
	print "Total correct-pair reads:\t$all_correct\n";
	print "Total mismatched-pair reads:\t$all_mismatched\n";
	print "Total PE1-only (expected):\t$all_pe1only_exp\n";
	print "Total PE1-only (other):\t$all_pe1only_oth\n";
	print "Total PE2-only (expected):\t$all_pe2only_exp\n";
	print "Total PE2-only (other):\t$all_pe2only_oth\n";
}
print "Total unidentified:\t$all_unidentified\n\n";
print "Grand total:\t$grand_total_reads\n\n\n";

print $log_fh "*" x 40 . "\n\n";
printf $log_fh "Total 3LTR read pairs:\t%s\n", exists $LTRs{'3LTR'} ? $all_3LTR : "N/A";
printf $log_fh "Total 5LTR read pairs:\t%s\n", exists $LTRs{'5LTR'} ? $all_5LTR : "N/A";
if ($use_PE2_inline) {
	print $log_fh "Total correct-pair reads:\t$all_correct\n";
	print $log_fh "Total mismatched-pair reads:\t$all_mismatched\n";
	print $log_fh "Total PE1-only (expected):\t$all_pe1only_exp\n";
	print $log_fh "Total PE1-only (other):\t$all_pe1only_oth\n";
	print $log_fh "Total PE2-only (expected):\t$all_pe2only_exp\n";
	print $log_fh "Total PE2-only (other):\t$all_pe2only_oth\n";
	# Per-combination breakdown of mismatched pairs.
	if ($all_mismatched > 0) {
		print $log_fh "\nMismatched-pair breakdown (likely index hopping):\n";
		print     "\nMismatched-pair breakdown (likely index hopping):\n";
		foreach my $LTR (sort keys %mismatched_totals) {
			foreach my $p1 (sort keys %{$mismatched_totals{$LTR}}) {
				foreach my $p2 (sort keys %{$mismatched_totals{$LTR}{$p1}}) {
					my $n = $mismatched_totals{$LTR}{$p1}{$p2};
					my $p1_name = $pe1_seq2name{$p1};
					my $p2_name = $pe2_seq2name{$p2};
					print     "\t$LTR\tPE1_${p1_name}_PE2_${p2_name}\t$n\n";
					print $log_fh "\t$LTR\tPE1_${p1_name}_PE2_${p2_name}\t$n\n";
				}
			}
		}
		print "\n";
		print $log_fh "\n";
	}
}
print $log_fh "Total unidentified pairs:\t$all_unidentified\n\n";
print $log_fh "Grand total:\t$grand_total_reads\n\n";
print "Compressing...\n\n\n";

foreach my $LTR (@{$settings{LTR}}) {
	my @files = (
		<$LTR/PE1*.txt>,
		<$LTR/PE2*.txt>,
		<$LTR/mismatched_reads/PE1*.txt>,
		<$LTR/single_index/PE1*.txt>,
		<$LTR/single_index/PE2*.txt>,
	);
	foreach my $f (@files) {
		if (! -s $f) {
			unlink $f or warn "Could not unlink empty file $f: $!\n";
			next;
		}
		if (-e $f) {
			my $rc = system("gzip", $f);
			if ($rc != 0) {
				my $signal = $rc & 127;
				my $exit   = $rc >> 8;
				die "gzip failed on $f (signal=$signal, exit=$exit): $!\n";
			}
		}
	}
}

# Print the time that the analysis completed to the log.
my $end_time = localtime;
print $log_fh "Analysis completed at $end_time\n\n";
print $log_fh "#" x 90 . "\n\n\n";
close $log_fh;
exit;

# This subroutine extracts the data from the fastq.gz files, identifies the LTR region 
# based on the sequences in the settings file, and writes reads to .txt files in
# subdirectories under <LTR>/. The data disappears after each pass through the subroutine.
#
# When $use_PE2_inline is enabled, reads are sorted into the following categories:
#   - PE1+PE2 both expected, valid pair  -> <LTR>/PE1_<exp>_PE2_<exp>.txt          (correct)
#   - Any combo where at least one side is "other" (foreign sample, likely hopping)
#                                        -> <LTR>/mismatched_reads/PE1_<X>_PE2_<Y>.txt
#   - Only PE1 matches (expected)        -> <LTR>/single_index/PE1_expected.txt
#   - Only PE1 matches (other)           -> <LTR>/single_index/PE1_other.txt
#   - Only PE2 matches (expected)        -> <LTR>/single_index/PE2_expected.txt
#   - Only PE2 matches (other)           -> <LTR>/single_index/PE2_other.txt
# Reads matching no known barcode (or failing the LTR/linker checks) are unidentified.
# When disabled, the script behaves as in the legacy version (PE1 match only,
# files written directly under <LTR>/).
sub extract_and_print_data {
	# Initialize the passed-in variables.
	my($pairs, $filenametemplate,
	   $file_correct, $file_mismatched,
	   $file_pe1only_exp, $file_pe1only_oth,
	   $file_pe2only_exp, $file_pe2only_oth,
	   $pe1_seq2name, $pe2_seq2name,
	   $pe1_is_expected, $pe2_is_expected,
	   $expected_pair,
	   $analyses, $match_LTR_seq, $log_fh,
	   $use_PE2_inline, $UMI_offset, $UMI_length) = @_;
	
	# Extract the contents of the fastq.gz files to temporary locations for processing.

	# Create the names of the files to extract using the file name template and the current $pairs.
	my ($R1, $R2);

	if (-e "${filenametemplate}R1.fastq") {
		$R1 = "${filenametemplate}R1.fastq";
		$R2 = "${filenametemplate}R2.fastq";

	} elsif (-e "${filenametemplate}R1.fastq.gz") {
		$R1 = "${filenametemplate}R1.fastq.gz";
		$R2 = "${filenametemplate}R2.fastq.gz";

	} elsif (-e "${filenametemplate}R1_${pairs}.fastq.gz") {
		$R1 = "${filenametemplate}R1_${pairs}.fastq.gz";
		$R2 = "${filenametemplate}R2_${pairs}.fastq.gz";

	} else {
		die "\nCan't find matching FASTQ files for template ${filenametemplate}\n";
	}
	
	# Print the file names to the monitor and the log.
	print "$R1\n$R2\n\n";
	print $log_fh "$R1\n$R2\n\n";

	# Open file handles for printing the sequence data.
	my ($INread1, $INread2);

	if ($R1 =~ /\.gz$/) {
		open($INread1, "gzip -dc $R1 |") or die "Unable to extract the read 1 sequences: $!";
		open($INread2, "gzip -dc $R2 |") or die "Unable to extract the read 2 sequences: $!"
	} else {
		open($INread1, "<", $R1) or die "Unable to extract the read 1 sequences: $!";
		open($INread2, "<", $R2) or die "Unable to extract the read 2 sequences: $!"
	}
	
	my $totalreads = my $_3LTR_all_reads = my $_5LTR_all_reads  = my $unidentified  = my $correct_count  = my $mismatched_count  = 0;my $pe1only_exp_count = my $pe1only_oth_count = my $pe2only_exp_count = my $pe2only_oth_count = 0;

	# Per-file counters for the log: $file_counts{$LTR}{$relative_path} = N
	my %file_counts;

	# Per-pass open output filehandles, lazily opened on first write:
	#   $fh_cache{$LTR}{$relative_path} = open filehandle (append mode)
	# This streams reads to disk as they arrive instead of buffering all of them in memory.
	my %fh_cache;

	# Helper closure that returns (and opens-on-first-use) the filehandle for a given LTR + rel_path.
	my $get_fh = sub {
		my ($LTR, $rel_path) = @_;
		return $fh_cache{$LTR}{$rel_path} if exists $fh_cache{$LTR}{$rel_path};
		my $file2print = "$LTR/$rel_path";
		open my $fh, ">>", $file2print
			or die "Can't open $file2print for printing: $!\n";
		$fh_cache{$LTR}{$rel_path} = $fh;
		return $fh;
	};

	# Per-combination mismatched-pair counter for this pass:
	my %mm_breakdown;
	
	
	# Begin unzipping and debarcoding the data
	SCAN: while (my $lineread1=<$INread1>) { # Begin with read 1
		chomp($lineread1);
	
		if($lineread1=~/^@(.*)( 1.*)/) { # Lines beginning with @ indicate sequence name 
			$totalreads++; 
			my $read1_id=$1; # Only want the part of the ID that's identical to read 2
			my $read1_seq=<$INread1>; # The next line is the sequence
			chomp($read1_seq);
			<$INread1>;  # Skip the + sign in between seq and quality
			my $read1_qual=<$INread1>; # Pull in the quality line
			chomp($read1_qual);
		
			# Repeat the same steps for read 2
			my $read2_id=<$INread2>;
			chomp($read2_id);
			$read2_id =~s/\@(.*) .*/$1/;
			my $read2_seq=<$INread2>;
			chomp($read2_seq);
			<$INread2>;
			my $read2_qual=<$INread2>;
			chomp($read2_qual);
			
			# Check to see if the sequence IDs match. Exit if they don't.
			if($read1_id ne $read2_id){
				print "Read 1 sequence ID does not match read 2:\n$read1_id\t$read2_id\n";
				print $log_fh "Read 1 sequence ID does not match read 2:\n$read1_id\t$read2_id\n";
				print "The read 1 and 2 sequence IDs don't match. There might be something wrong with your data.\n";
				print "All sequence ID must match to complete the analysis.";
				exit;
			}
			
			# Match the sequence to 3LTR, 5LTR, or both, depending on user input.
			my ($R1_junction, $LTR) = &$match_LTR_seq($read1_seq);

			# If neither nested primer was detected, skip to the next sequence.
			unless (defined $LTR) {
				$unidentified++;
				next SCAN;
			}
			
			# --- PE1 inline barcode (start of Read 1, after the 6nt dogtag) ---
			my $read1_bc = substr($read1_seq, 6, 8);
			my $pe1_match = exists $pe1_seq2name->{$read1_bc};

			# --- PE2 inline barcode (start of Read 2, only when enabled) ---
			my $pe2_match = 0;
			my $read2_bc;
			if ($use_PE2_inline) {
				$read2_bc = substr($read2_seq, 0, 8);
				$pe2_match = exists $pe2_seq2name->{$read2_bc};
			}

			# Reads with no matching barcode are unidentified.
			if ($use_PE2_inline) {
				if (!$pe1_match and !$pe2_match) { $unidentified++; next SCAN; }
			} else {
				if (!$pe1_match) { $unidentified++; next SCAN; }
			}

			# Only compare the LTR and barcode sequences to the first 75 bases of the reads.
			my $read1_first75=substr $read1_seq,0,75; 
			my $read2_first75= substr $read2_seq,0,75;
			
			# Store the data for later printing if the junction and linker sequences.
			if	($read1_first75 =~ m/$R1_junction/ and $read2_first75 =~ m/TCCGCTTAGAGGACT/) {
				# Keep track of the number of 3LTR and 5LTR reads.
				if ($LTR eq "3LTR") {
					$_3LTR_all_reads++;
				} else {
					$_5LTR_all_reads++;
				}
				# Split read 1 based on the HIV/host junction sequence.
				my ($preLTRend, $postLTRend) = split (/$R1_junction/, $read1_seq);
				my $LTR_junction= "$preLTRend$R1_junction\t$postLTRend";
				# Get the molecular identifier (MID) FIRST, from the unmodified read 2:
				# bases 0-9 (legacy) or 8-17 (PE2 inline mode).
				my $MID = substr $read2_seq, $UMI_offset, $UMI_length;
				# Now insert a tab at the end of the linker sequence in read 2.
				my ($preLinker, $postLinker) = split (/TCCGCTTAGAGGACT/, $read2_seq);
				$read2_seq = "$preLinker"."TCCGCTTAGAGGACT\t$postLinker";

				# Choose destination file (relative path under <LTR>/) based on which barcodes matched
				# and whether each is "expected" or "other".
				my $rel_path;
				if ($use_PE2_inline) {
					if ($pe1_match and $pe2_match) {
						my $pe1_exp = exists $pe1_is_expected->{$read1_bc};
						my $pe2_exp = exists $pe2_is_expected->{$read2_bc};
						if ($pe1_exp and $pe2_exp
							and defined $expected_pair->{$read1_bc}
							and $expected_pair->{$read1_bc} eq $read2_bc) {
							# Both expected and form the correct pair.
							$rel_path = $file_correct->{$LTR}{$read1_bc};
							$correct_count++;
						} else {
							# At least one side is "other" (or expected/expected but wrong combo).
							my $fname = $file_mismatched->{$LTR}{$read1_bc}{$read2_bc};
							$rel_path = "mismatched_reads/$fname" if defined $fname;
							$mismatched_count++;
							$mm_breakdown{$LTR}{$read1_bc}{$read2_bc}++;
						}
					} elsif ($pe1_match) {
						if (exists $pe1_is_expected->{$read1_bc}) {
							$rel_path = "single_index/" . $file_pe1only_exp->{$LTR};
							$pe1only_exp_count++;
						} else {
							$rel_path = "single_index/" . $file_pe1only_oth->{$LTR};
							$pe1only_oth_count++;
						}
					} else { # pe2_match only
						if (exists $pe2_is_expected->{$read2_bc}) {
							$rel_path = "single_index/" . $file_pe2only_exp->{$LTR};
							$pe2only_exp_count++;
						} else {
							$rel_path = "single_index/" . $file_pe2only_oth->{$LTR};
							$pe2only_oth_count++;
						}
					}
				} else {
					# Legacy mode: directly under <LTR>/, named PE1_<name>_PE2_<single name>.txt.
					$rel_path = $file_correct->{$LTR}{$read1_bc};
				}

				unless (defined $rel_path) { $unidentified++; next SCAN; }

				my $out_fh = $get_fh->($LTR, $rel_path);
				print $out_fh "$read1_id#$MID\t1\t$LTR_junction\t$read1_qual\t$read2_id#$MID\t2\t$read2_seq\t$read2_qual\t$MID\n"
					or die "Write to $LTR/$rel_path failed: $!\n";
				$file_counts{$LTR}{$rel_path}++;
			}else {$unidentified++;}
		}
	}		
	close $INread1;
	close $INread2; 
	
	# Print out the sequence data from this pass through the subroutine to files and print the 
	# read counts to the terminal and log file. 
	foreach my $LTR (sort keys %fh_cache) {
		foreach my $rel_path (sort keys %{$fh_cache{$LTR}}) {
			close $fh_cache{$LTR}{$rel_path}
				or die "Close failed on $LTR/$rel_path: $!\n";
		}
	}
	
	# Print per-file counts to the terminal and log file. Data is already written.
	foreach my $LTR (sort keys %file_counts) {
		print "$LTR\n";
		print $log_fh "$LTR\n";
		
		foreach my $rel_path (sort keys %{$file_counts{$LTR}}) {
			my $counts = $file_counts{$LTR}{$rel_path} // 0;
			my $label;
			if ($rel_path =~ m{^mismatched_reads/}) {
				$label = "mismatched-pair reads";
			} elsif ($rel_path =~ m{^single_index/PE1_expected\.txt$}) {
				$label = "PE1-only reads (expected)";
			} elsif ($rel_path =~ m{^single_index/PE1_other\.txt$}) {
				$label = "PE1-only reads (other)";
			} elsif ($rel_path =~ m{^single_index/PE2_expected\.txt$}) {
				$label = "PE2-only reads (expected)";
			} elsif ($rel_path =~ m{^single_index/PE2_other\.txt$}) {
				$label = "PE2-only reads (other)";
			} elsif ($rel_path =~ /^PE1_.+_PE2_.+\.txt$/) {
				if ($use_PE2_inline) {
					$label = "correct-pair reads";
				} else {
					$label = "read pairs";
				}
			} else {
				$label = "read pairs";
			}
			
			print     "\t$rel_path\n\t\t$LTR $label:\t$counts\n\n";
			print $log_fh "\t$rel_path\n\t\t$LTR $label:\t$counts\n\n";
		}
	}

	# Print the unidentified read counts to the terminal and the log file.
	print "Unidentified pairs:\t$unidentified\n\n";
	print $log_fh "Unidentified pairs:\t$unidentified\n\n";
	return ($totalreads, $_3LTR_all_reads, $_5LTR_all_reads, $unidentified,
        $correct_count, $mismatched_count,
        $pe1only_exp_count, $pe1only_oth_count,
        $pe2only_exp_count, $pe2only_oth_count,
        \%mm_breakdown);	
}

###################################################################################################

# This subroutine returns a subroutine reference that will be used to analyze the LTR data
sub choose_LTR_sub {
	my %IUPAC = (  # This hash replaces the ambiguous bases in the primer with a regex
	"R"	=> "[AG]", # character class for downstream matching and substitutions.
	"Y"	=> "[CT]",
	"S"	=> "[GC]",
	"W"	=> "[AT]",
	"K"	=> "[GT]",
	"M"	=> "[AC]",
	"B"	=> "[CGT]",
	"D"	=> "[AGT]",
	"H"	=> "[ACT]",
	"V"	=> "[ACG]",
	"N" => "[ATGC]",
	);
	
	my($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction, @analyses)=@_; # Initialize the passed in variables.
	$_ = uc $_ foreach ($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction); # Make all base symbols uppercase if they're not already
	my $_3LTR_primer_length = length($_3LTR_primer) if defined $_3LTR_primer and $_3LTR_primer ne '';
	my $_5LTR_primer_length = length($_5LTR_primer) if defined $_5LTR_primer and $_5LTR_primer ne '';
	# Replace the ambiguous bases in the primers and junctions with character classes.
	foreach my $seq ($_3LTR_primer, $_3LTR_junction, $_5LTR_primer, $_5LTR_junction) {
		next unless defined $seq and $seq ne '';
		$seq =~ s/([RYSWKMBDHVN])/$IUPAC{$1}/g;
	}
	
	# Decide which analysis subroutine to use.  
	my $num_choices = @analyses; # Store the number of items in @analyses as a scalar variable.
	foreach (@analyses) {croak "Invalid LTR option selected. Check your settings file.\n" unless m/^[35]LTR$/;} # Croak if anything other than 3LTR or 5LTR is in %analyses.
	croak "More than 2 options were entered for LTR analysis. Check your settings file.\n" if $num_choices > 2;  # Croak if more than 2 items were entered for LTR analyses.
	croak "You entered $analyses[0] twice. Check your settings file." if $num_choices == 2 and $analyses[0] eq $analyses[1]; # Croak if the same LTR was entered twice.
	if ($num_choices == 2) {
		return sub {
			my ($read1_seq) = @_;
			my $read1_nestedPrimer;
			if ($read1_nestedPrimer=substr($read1_seq,14,$_3LTR_primer_length) and $read1_nestedPrimer =~ /$_3LTR_primer/) {
				my ($R1_junction) = "$_3LTR_junction";
				return($R1_junction, "3LTR");
			
			}elsif ($read1_nestedPrimer=substr($read1_seq,14,$_5LTR_primer_length) and $read1_nestedPrimer =~ /$_5LTR_primer/) {
				my ($R1_junction) = "$_5LTR_junction";
				return($R1_junction, "5LTR");
			
			}else {
				return;
			}
		};
		
	# Use the same steps above but only check for a match to the 3LTR primer.
	}elsif ($analyses[0] eq "3LTR") {
		return sub {
			my ($read1_seq) = @_;
			my $read1_nestedPrimer=substr($read1_seq,14,$_3LTR_primer_length);
			return unless $read1_nestedPrimer =~ /$_3LTR_primer/;
			my ($R1_junction) = "$_3LTR_junction";
			return($R1_junction, "3LTR");
		};
		
	# Use the same steps above but only check for a match to the 5LTR primer.	
	}elsif ($analyses[0] eq "5LTR") {
		return sub {
			my ($read1_seq) = @_;
			my $read1_nestedPrimer=substr($read1_seq,14,$_5LTR_primer_length);
			return unless $read1_nestedPrimer =~ /$_5LTR_primer/;
			my ($R1_junction) = "$_5LTR_junction";
			return($R1_junction, "5LTR");
		};
	}
}