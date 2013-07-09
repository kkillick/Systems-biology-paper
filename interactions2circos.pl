use warnings;
use strict;
use File::Basename;
use Getopt::Long;

my $USAGE = "

# Karsten Hokamp (kahokamp\@tcd.ie), Trinity College Dublin, 2013
# generate input files for Circos from InnateDB interactions
# and lists of differentially expressed genes

USAGE: $0 [options] interaction_file gene_lists

Options:
  -ref    name of gene used as centre piont
  -prefix prefix for output files
  -fdr    FDR cut-off
  -help   this text

interaction_file: generated through InnateDB
example: 
Group ID        Interaction ID(s)       Interaction     Interactors(Symbol|InnateDB ID|Ensembl ID|Entrez Gene ID|Type)  Interactor Types        Interactor Species      Interaction Level       Interaction Type        Source Database ID(s)   Supporting Publications PMID    Evidences(IDB-ID|Experiment type|Tissue/cell type|PMID) Human Orthologs(Symbol|IDBG-ID|Ensembl) Mouse Orthologs(Symbol|IDBG-ID|Ensembl) Bovine Orthologs(Ensembl)
55233   IDB-224010      Transcription factor EGR1 binds MIR455 gene     EGR1|IDBG-47191|ENSG00000120738|1958|protein :: MIR455|IDBG-126451|ENSG00000207726|null|dna      protein - DNA  Homo sapiens    direct interaction      physical association    IDB-224010      1       20811575        IDB-224010|chromatin immunoprecipitation assay|chronic myeloid leukemia cell line|20811575  ;           Egr1|IDBG-136414|ENSMUSG00000038418 :: [no MIR455 ortholog found]       EGR1|IDBG-645986|ENSBTAG00000010069 :: [no MIR455 ortholog found]
...

gene_lists: 
example:
ID      GeneSymbol      Gene    mismatch        ENS human       ENS cow GeneDesc        Coef.(BOV_24Hr - BCG_24Hr)      t.(BOV_24Hr - BCG_24Hr) p.value.(BOV_24Hr - BCG_24Hr)   p.value.adj.(BOV_24Hr - BCG_24Hr)
Bt.15822.1.S1_at        FGR     FGR             ENSG00000000938 ENSBTAG00000011784      Gardner-Rasheed feline sarcoma viral (v-fgr) oncogene homolog   0.490641027     3.625181919     0.000511322     0.003221085
(GeneSymbol, Coef and p.value.adj are extracted)

";

my $help = '';

my $ref = '';
my $prefix = 'circos';
my $fdr_threshold = 0.01;

&GetOptions(
    'fdr_threshold=f' => \$fdr_threshold,
    'ref=s' => \$ref,
    'prefix=s' => \$prefix,
);


# Step 1: read in InnateDB interactions

my $innatedb = shift;
open (IN, $innatedb)
    or die;
my %innatedb = ();
my $header = <IN>;
my %headers = &get_headers($header);

while (<IN>) {
    chomp;
    my @h = split /\t/, $_;
    my $ids = '';
    if (defined $headers{'Interactors(Symbol|InnateDB ID|Ensembl ID|Entrez Gene ID|Type)'}) {
	$ids = $h[$headers{'Interactors(Symbol|InnateDB ID|Ensembl ID|Entrez Gene ID|Type)'}];
    } elsif (defined $headers{'Interactors (Symbol|InnateDB ID|Ensembl ID|Entrez Gene ID|Type)'}) {
	$ids = $h[$headers{'Interactors (Symbol|InnateDB ID|Ensembl ID|Entrez Gene ID|Type)'}];
    }

    foreach my $id (split /\:\:/, $ids) { #/) {
	next unless ($id =~ /^\s*(\w.*?)\|.+\|(ENSG\d+)/);
	my $symbol = $1;
	my $ensg = $2;
	$innatedb{$symbol}{$ensg}++;
    }
}
close IN;

my $num = scalar keys %innatedb;

print STDERR "Read $num interactors from $innatedb\n";


# Step 2: read in differentially expressed genes:
my %in = ();
foreach my $file (@ARGV) {
    unless ($file =~ /(\d+)h.?_paired/) {
	die "Wrong file: $file\n";
    }
    my $hour = $1;
    if (length($hour) < 2) {
	$hour = '0'.$hour;
    }

    open (IN, $file) 
	or die;
    my $header = <IN>;

#    0  'ID'
#   1  'GeneSymbol'
#   2  'Gene'
#   3  'mismatch'
#   4  'ENS human'
#   5  'ENS cow'
#   6  'GeneDesc'
#   7  'Coef.(BOV_24Hr - BCG_24Hr)'
#   8  't.(BOV_24Hr - BCG_24Hr)'
#   9  'p.value.(BOV_24Hr - BCG_24Hr)'
#   10  "p.value.adj.(BOV_24Hr - BCG_24Hr)\cM\cJ"

    while (<IN>) {
	chomp;
	s/\cM//;
	my @h = split /\t/, $_;
	push @{$in{$hour}{$h[1]}{$h[10]}}, $h[7];
    }
    close IN;
}


# Step 3: merge values from multiple instances of the same gene:

my %in2 = ();
my %genes = ();

foreach my $hour (sort keys %in) {
    foreach my $gene (sort keys %{$in{$hour}}) {
	my @fdr = (sort { $a <=> $b } keys %{$in{$hour}{$gene}});
	my @fc = @{$in{$hour}{$gene}{$fdr[0]}};
	my $fc = &median(@fc);
	$in2{$hour}{$gene}{'fdr'} = $fdr[0];
	$in2{$hour}{$gene}{'fc'} = $fc;
	if ($fdr[0] > $fdr_threshold) {
	    next;
	}
	$genes{$gene}++;
    }
}
%in = %in2;


# Step 4: print labels file:

open (OUT, ">$prefix.labels")
    or die;
my $start = -1;
foreach my $gene (sort keys %genes) {
    next if ($gene eq $ref);
    next unless (defined $innatedb{$gene});

    $start++;
    my $end = $start + 1;;
    print OUT "chr1 $start $end $gene\n";
}
close OUT;


# Step 5: write data files

foreach my $hour (sort keys %in) {
    open (OUT, ">$prefix.$hour.txt")
	or die;
    my $start = -1;
    foreach my $gene (sort keys %genes) {
	next if ($gene eq $ref);
	next unless(defined $innatedb{$gene});

	my $val = 0;
	unless (defined $in{$hour}{$gene}) {
	    warn "No value for $hour, $gene\n";
	} else {
	    $val = $in{$hour}{$gene}{'fc'};
	}
	$start++;
	my $end = $start + 1;
	print OUT "chr1 $start $end $val\n";
    }
    close OUT;
}

open (OUT, ">$prefix.combined.txt")
    or die;
select OUT;
print "Gene";
foreach my $hour (sort keys %in) {
    print "\t$hour";
}
print "\n";

foreach my $gene (sort keys %genes) {
    next unless(defined $innatedb{$gene});

    my @out = ($gene);

    foreach my $hour (sort keys %in) {
	
	my $val = '';
        if (defined $in{$hour}{$gene}) {
            $val = $in{$hour}{$gene}{'fc'};
        }
	push @out, $val;
    }

    print "".(join "\t", @out)."\n";
}
close OUT;


# Subroutines:
    
sub median {

    unless (@_) {
	return '';
    }

    if (@_ == 1) {
	return $_[0];
    }

    my @order = sort { $a <=> $b } @_;
    my $md = '';

    if ($#order % 2 == 0){
        $md = $order[$#order / 2];
    } else {
	my $index1 = int($#order / 2);
	my $index2 = int($#order / 2) + 1;
	my $low = $order[$index1];
	my $high = $order[$index2];
	$md = $low + ($high - $low)/2;
    }
    return $md;
}

sub get_headers {
    my %out = ();
    my $header = shift;
    chomp $header;
    my $count = 0;
    my @tmp = split /\t/, $header;
    foreach (@tmp) {
	$out{$_} = $count;
	$count++;
    }
    return %out;
}
