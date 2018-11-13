use strict;
print "$0";
die "usage: perl program.pl <infile.anno>  <output.tsv>" unless @ARGV==2;

my $file1 = shift;
my $out = shift;


trans2tsv($file1,$out);


sub trans2tsv{
	my $f = shift;
	my $of = shift;
	open IN  ,"<$f" || die $!;
	if(-e $of){
		open OUT,">>$of" || die $!;
	}
	else{
		open OUT,">$of" || die $!;
		print OUT "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\tvariant_case\tvariant_freq\tgenotype\n";
	}
	<IN>;
	while(<IN>){
		s/[\r\n]+//;
		my @cells = split /\t/;
		#next if $cells[5] ne 'exonic' and $cells[5] ne 'splicing';
		my $mid = getmid($cells[5],$cells[9],$cells[7],"$cells[0]\_$cells[1]\_$cells[2]\_$cells[3]\_$cells[4]");
		my ($rc,$vc,$vf);
		my $nc = 2;
		my $minor_cn = 0;
		my $major_cn = 2;
		my $variant_case = $cells[-1];
		my $genotype = 'unknown';
		($rc,$vc,$vf) = getinfo($cells[-1]);
		next if ($rc +$vc) < 50;
		print OUT "$mid\t$rc\t$vc\t$nc\t$minor_cn\t$major_cn\t$variant_case\t$vf\t$genotype\n";
	}
	close IN;
	close OUT;
}

#0/1:245,13:0.052:9:4:.:8288,434:139:106

sub getinfo{
	my $info = shift;
	my $trc=1;
	my $tvc=1;
	my $tvf=1;
	if($info =~/[\d]\/[\d]:(\d+),(\d+):/){
		$trc = $1;
		$tvc = $2;
		$tvf = $tvc/($trc+$tvc);
	}
	return ($trc,$tvc,$tvf);
}


sub getmid{
	my $type = shift;
	my $name1 = shift;
	my $name2 = shift;
	my $name3 = shift;
	if($type eq 'exonic'){
		return $name1;
	}
	elsif($type eq 'splicing'){
		return $name2;
	}
	else{
		return $name3;
	}
}
