use strict;

die "usage: perl prog <table.old_style> <HBXXXX-5-Y.txy> <output>" unless @ARGV ==3;


my $table = shift;
my $driver = shift;
my $out = shift;


#NM\_002520(NPM1):c.859\_860insTCTG(p.L287fs)	Exon 11	25.39\%	allgene	xyb1NPM001


my %hash;
open IN,$driver || die $!;
while(<IN>){
	s/[\r\n]+//;
	my $gene;
	my $pdot;
	if(/\((\w+)\).+?\((p.+)\)/){
		$gene = $1;
		$pdot = $2;
	}
	elsif(/\((\w+)\).+?[+-]/){
		$gene = $1;
		$pdot = 'splicing';
	}
	$hash{$gene}{$pdot} = 1 if $gene ne '';
}
close IN;

open IN,$table || die $!;
open OUT,">$out" || die $!;
my $title = <IN>;
$title =~ s/cluster_id/cluster/;
print OUT "Driver\tid\t$title";
while(<IN>){
	s/[\r\n]+//;
	my @cells = split /\t/;
	my $pass =1;
	foreach my $gene(keys %hash){
#	    print "gene:$gene\n";
		if($cells[0] =~/$gene/){
			foreach my $pd(keys %{$hash{$gene}}){
				if($cells[0]=~/$pd/){
					print OUT "TRUE\t$gene\_$pd\t$_\n";
					$pass=0;
				}
				elsif($pd eq 'splicing' and $cells[0]=~/[+-]/){
					print OUT "TRUE\t$gene\_$pd\t$_\n";
					$pass =0;
				}
			}
		}
	}
	if($pass == 1){
		print OUT "FALSE\t-\t$_\n";
	}
}
close OUT;
