#!/usr/bin/perl -w
$file =$ARGV[0];
$file2=$ARGV[1];
$file3=$ARGV[2];
open(FH,$file);
open(FO1,">".$file2);
open(FO2,">".$file3);
@aa= ('A' ,'C' ,'D' ,'E' ,'F' ,'G' ,'H' ,'I' ,'K' ,'L' ,'M' ,'N' ,'P' ,'Q' ,'R' ,'S' ,'T' ,'V' ,'W' ,'Y');
while ($line=<FH>){ 
	chomp $line;$line=~s/\s//g;$line=uc($line);
	@peptide=split('',$line);
	for($i=0;$i<@peptide;$i++){
		for($j=0;$j<20;$j++){
		        undef @new;
		        undef @acc;
			@new=@peptide;
		        @acc=@peptide;
			if($new[$i] ne $aa[$j]){
			    #$acc[$i]='<font color =red>'.$aa[$j].'</font>';
			    $new[$i]=$aa[$j];			
			    $newline=join ('',@new);
			    #$newacc=join ('',@acc);
			    $pos=$i+1;
			    print FO1 $newline."\t".$pos."\n";
			    #print FO2 $newacc."\n";
			}
			if($acc[$i] ne $aa[$j]){
			    $acc[$i]='<font color =red>'.$aa[$j].'</font>';
			   $newacc=join ('',@acc); 
			  print FO2 $newacc."\n";
			} 
		}
	}
}close FO1;
close FO2;
