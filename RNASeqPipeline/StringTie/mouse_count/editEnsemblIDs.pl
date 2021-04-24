$gtf = $ARGV[0];

open(IN, $gtf) || die $!;
open(OUT, ">temp.gtf") || die $!;

while(<IN>){
	if(/(ENSMUSG\d+)\.\d+/){
		$rep = $1;
		s/ENSMUSG\d+\.\d+/$rep/;
		}
	if(/(ENST\d+)\.\d+/){
		$rep = $1;
		s/ENST\d+\.\d+/$rep/;
		}
	print OUT;
	}

close IN;
close OUT;

`mv temp.gtf $gtf`;
 
