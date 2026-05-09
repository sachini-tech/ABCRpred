#!/usr/bin/perl
$file=$ARGV[0];			##----sequence file can be multiple sequneces no fasta header----------##

$file_out=$ARGV[1];   		##------------------output file----------------##
$file_out2=$ARGV[2];
$properties="/usr1/webserver/cgidocs/raghava/il4pred/prop";
open(FH,$file);
open(FP, ">$file_out");
open(FP2, ">$file_out2");
$count=1000;
while ($seq=<FH>){$count++;
chomp($seq); 
	#print FP "$count\t";
	#print FP "$seq\t";
		  %hash= stop($seq,"$properties/hydrophobicity.txt");  
		  printf (FP "<td align = center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash= stop($seq,"$properties/steric.txt");
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash= stop($seq,"$properties/sidebulk.txt");
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash= stop($seq,"$properties/hydrpathy.txt");
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash= stop($seq,"$properties/amphipathicity.txt");
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash1= stop($seq,"$properties/hydrophilicity.txt");
		  #printf (FP2 "%3.2f\t",$avvg);
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash= stop($seq,"$properties/nethydrogen.txt");
		  printf (FP "<td align=center>%3.2f</td>",$avvg);
		  printf (FP2 "%3.2f\t",$avvg);
		  %hash11= stop($seq,"$properties/charge.txt");
		  printf (FP "<td align=center>%3.2f</td>",$sum);
		  printf (FP2 "%3.2f\t",$sum);
		  %hash12= pi($seq);
		  printf (FP "<td align=center>%3.2f</td>",$ph);
		  printf (FP2 "%3.2f\t",$ph);			
@ab=split('',uc($seq)); $n=18*(@ab-1);
                $val=0;$sum=0;
                        foreach $res(@ab){
                                open (FH1,"$properties/mol_wt.txt");@mat=();
                                while ($line1=<FH1>){
                                        chomp $line1;
                                        if ($line1 =~ /$res/) {
                                                @mat=split"\t",$line1;
                                                $val=$mat[1];
                                            }
                                }
                                $sum +=$val;
                        }
                        $net=($sum-$n);
                printf (FP "<td align=center>%3.2f</td></tr>\n",$net);
		  printf (FP2 "%3.2f\n",$net);
}



##-----------subroutine starts-------#######

sub stop{
    ($seq,$path)=@_;
    $seq=uc($seq);$len=length($seq);
    @seq=split("",$seq);
    open(PRO,$path);
    @prop=<PRO>;
    close PRO;
    foreach $line(@prop){
        chomp($line);
	@split=split("\t",$line);
        $property{$split[0]}=$split[1];
    }
    $sum=0;$avg=0;
    $aa = "ACDEFGHIKLMNPQRSTVWYX";
    @aa=split("",$aa);
    foreach(  @aa){$sum_p{$_}=0;}
    foreach $res(@seq){
        $val=$property{$res};
        $sum_p{$res}=$sum_p{$res}+$val;
        $sum += $val;
    }$avvg=($sum/$len);#print "$avvg\n";
    return %sum_p;
}

sub pi{
	($seq,$path)=@_;
	$seq=uc($seq);
	@pep=split('',$seq);
        #print"$pep[0]\n";
        $asp=$glu=$cys=$tyr=$his=$lys=$arg=0;
                for($i=0;$i<=@pep;$i++){
                if($pep[$i] eq 'D'){$asp++;}
                elsif($pep[$i] eq 'E'){$glu++;}
                elsif($pep[$i] eq 'C'){$cys++;}
                elsif($pep[$i] eq 'Y'){$tyr++;}
                elsif($pep[$i] eq 'H'){$his++;}
                elsif($pep[$i] eq 'K'){$lys++;}
                elsif($pep[$i] eq 'R'){$arg++;}
                }
        #print"D\t$asp\nE\t$glu\nC\t$cys\nY\t$tyr\nH\t$his\nK\t$lys\nR\t$arg\n";
        $ph=0;
        do{$ph += 0.01;#print"$ph\n";
        $qn1=$qn2=$qn3=$qn4=$qn5=$qp1=$qp2=$qp3=$qp4=0;
        $qn1=-1/(1+10**(3.55-$ph));             #C-terminal charge
        $qn2=-$asp/(1+10**(4.05-$ph));           #D charge
        $qn3=-$glu/(1+10**(4.45-$ph));            #E charge
        $qn4=-$cys/(1+10**(9-$ph));            #C charge
        $qn5=-$tyr/(1+10**(10-$ph));        #Y charge
        $qp1=$his/(1+10**($ph-5.98));            #H charge
        $qp2=1/(1+10**($ph-8.2));                #NH2charge
        $qp3=$lys/(1+10**($ph-10));            #K charge
        $qp4=$arg/(1+10**($ph-12));            #R charge
        $nq=$qn1+$qn2+$qn3+$qn4+$qn5+$qp1+$qp2+$qp3+$qp4;
        }
        until($nq <= 0);
	return %ph;
}
