#!/usr/bin/perl -w
$file = $ARGV[0];	####-------input multipeptide file------####
$file2= $ARGV[1];	####-------output file------####
$id=$ARGV[2];		####-------first id value------####
open (FH, "$file");
open (FP, ">$file2");
while ($line=<FH>){
	chomp $line;
	$line=">$id ##".$line;
	print (FP "$line\n");$id++;
}
close FH;close FP;
exit;	
