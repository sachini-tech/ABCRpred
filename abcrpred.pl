#!/gpsr/local/bin/perl
use Getopt::Std;
getopts('i:o:t:');

$input=$opt_i;
$output=$opt_o;
$thresh=$opt_t;

if($opt_i eq '')
{
    print "USAGE: abcrpred.pl -i <fasta format sequences> -o <output file name> -t <threshold>\n\n";
    print "Example Command: ./abcrpred.pl -i /gpsr/examples/example_abcrpred.fasta -o out -t 0.53\n\n";
    print "-i\tSequence in FASTA format\n";
    print "-o\toutput file\n";
    print "-t\tthreshold\n";
    exit();
}
$prog_path="/gpsr/standalone/abcrpred/prog";
$ran=int(rand(10000));
$dir="/gpsr/standalone/abcrpred/temp/abcrpred_$ran";
$modeldir = "/gpsr/models/abcrpred";
`mkdir "$dir"`;
`/gpsr/software/anaconda3/bin/python3 $prog_path/read_fasta.py $input $dir/input_user_1`;
`grep -v '>' $dir/input_user_1 >$dir/input_sline`;
`grep '>' $dir/input_user_1 >$dir/header`;
`perl -pi -e 's/\r$//g' $dir/input_sline`;
if($method==1 or $method==""){
`/gpsr/software/anaconda3/bin/python3 $prog_path/abcrpred_feature_generation.py $dir/input_sline $dir/Out21 $dir`;
`/gpsr/software/anaconda3/bin/python3 $prog_path/top_10.py  $prog_path/top_33.csv $dir/Out21 $dir/top_33_features`;
`/gpsr/software/anaconda3/bin/python3 $prog_path/ml_run.py $dir/top_33_features $dir/Out3 $modeldir/RF_model`;
}
`/gpsr/software/anaconda3/bin/python3 $prog_path/determine.py $dir/Out3 $thresh $dir/Out333`;
`paste $dir/Out3 $dir/Out333 >$dir/Out66`;
`perl $prog_path/property_modified.pl $dir/input_sline $dir/Out6 $dir/Out44`;
`paste $dir/header $dir/Out66 $dir/Out44 > $dir/Out_sam`;
`printf '%s' "Seq_ID\tPrediction\tScore\tHydrophobicity\tSteric Hindrance\tSide bulk\tHydropathicity\tAmphipathicity\tNet Hydrogen\tCharge\tpI\tMol Wt.\n" | cat - $dir/Out_sam > $dir/RESULTS.csv`;

`cp $dir/RESULTS.csv $output`;

