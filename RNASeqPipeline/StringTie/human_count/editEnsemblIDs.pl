$gtf = $ARGV[0];

open(IN, $gtf) || die $!;
open(OUT, ">temp.gtf") || die $!;

while(<IN>){
	if(/(ENSG\d+)\.\d+/){
		$rep = $1;
		s/ENSG\d+\.\d+/$rep/;
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
 
