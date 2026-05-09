#!/usr/bin/perl -w
#perl code to create all possible mutants (with a user-defined window size) of a given peptide or list of peptides.
#usage:perl program.pl input_file output_file window_size 
$file =$ARGV[0];#input file--
$file2=$ARGV[1];#output file--
$sub=$ARGV[2];#window size--5 (for 5mer) or 10 for (10mer)
open(FH,$file);
open(FO,">".$file2);
while ($line=<FH>){
	chomp $line;$line=~s/\s//g;$line=uc($line);
	$len=length($line);
	$diff=$len-$sub;
	for ($t=0;$t<=$diff;$t++){
		$subline=substr($line,$t,$sub);
		print FO $subline."\n";
	}
}
