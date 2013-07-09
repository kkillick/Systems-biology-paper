use warnings;
use strict;
use Getopt::Long;

# list of hub/bottleneck genes to link connections to:
my $hubs = 'MYC,IKBKE,NFKB1';

my $USAGE = "

# Karsten Hokamp (kahokamp\@tcd.ie), Trinity College Dublin, 2013
# overlapping various degrees of genes connected to a set of 
# hub/bottlenecks with a list of differentially expressed genes

USAGE: $0 [options] gene_list_file cytoscape_file

Options:
  -hubs  list of gene names used as hubs, separated by comma
         default: '$hubs'
  -help  this text

gene_list_file: output from limma containing differentially expressed genes
example: 
line   A       GeneSymbol      GeneDesc        Coef.BOV_24Hr - BCG_24Hr        p.value.BOV_24Hr - BCG_24Hr     p.value.adj.BOV_24Hr - BCG_24Hr F       F.p.value       ID
1      6.95    SIGLEC1 \"sialic acid binding Ig-like lectin 1, sialoadhesin\"    3.137   0       0       43.32   0       Bt.6038.1.S1_at
...
(column  3 (GeneSymbol) will be extracted)

cytoscape_file: copy of all selected interactions (edge view) from Cytoscape, generated from within InnateDB
example:
http://www.innatedb.com:80/getInteractionCard.do?id=284316&group=true   physical association    A2M (physical association) IFIT5        A2M (physical association) IFIT5        A2M (physical association) IFIT5        1228.0  physical association    1
http://www.innatedb.com:80/getInteractionCard.do?id=284326&group=true   physical association    A2M (physical association) IL6  A2M (physical association) IL6  A2M (physical association) IL6  990.87263872    physical association    1
(interactions are extracted from column 4

";

my $help = '';

&GetOptions(
    'hubs=s' => \$hubs,
    'help|h' => \$help,
    );

if ($help) {
    print $USAGE;
    exit;
}

my %hubs = ();
foreach (split /, ?/, $hubs) { #/) {
    $hubs{$_}++;
}
my @hubs = sort keys %hubs;
my $hub_count = scalar @hubs;

print STDERR "Using the following $hub_count genes as hubs: ".(join ", ", @hubs)."\n";


# Step 1: read in file with list of differentially expressed genes:

my $list_file = shift;
open (IN, $list_file)
    or die;

# example lines from file with differentially expressed genes:
#line	A	GeneSymbol	GeneDesc	Coef.BOV_24Hr - BCG_24Hr	p.value.BOV_24Hr - BCG_24Hr	p.value.adj.BOV_24Hr - BCG_24Hr	F	F.p.value	ID
#1	6.95	SIGLEC1	"sialic acid binding Ig-like lectin 1, sialoadhesin"	3.137	0	0	43.32	0	Bt.6038.1.S1_at

my $head = <IN>;
my %genes = ();
my $i = 0;

# collect genes in hash, make them all upper case
while (<IN>) {
    $i++;
    my @h = split /\t/, $_;
    my $gene = uc $h[2];
    $genes{$gene}++;
}
close IN;

my $num = scalar keys %genes;

print STDERR "Extracted $num genes from $i entries in $list_file\n";


# Step 2: read all network connections from InnateDB

my %interactions = ();
my %connections = ();
my $cytoscape_file = shift;
my $all_connections = 0;
my $deg_connections = 0;
open (IN, $cytoscape_file)
    or die;
while (<IN>) {

# example lines from input file
#http://www.innatedb.com:80/getInteractionCard.do?id=284316&group=true   physical association    A2M (physical association) IFIT5        A2M (physical association) IFIT5        A2M (physical association) IFIT5        1228.0  physical association    1
#http://www.innatedb.com:80/getInteractionCard.do?id=284326&group=true   physical association    A2M (physical association) IL6  A2M (physical association) IL6  A2M (physical association) IL6  990.87263872    physical association    1

    my @h = split /\t/, $_; 
    $h[3] =~ s/ gene//ig;
    my ($interactor1, $interactor2) = split /\s+\(.+\)\s+/, $h[3]; 
    $interactor1 = uc $interactor1; 
    my @tmp = split / \:\: /, $interactor1; # /; 
    my @interactor1 = ();
    foreach my $node (@tmp) {
	foreach (split /;/, $node) { #/) {
	    push @interactor1, $_;
	}
    }

    $interactor2 = uc $interactor2; 
    @tmp = split / \:\: /, $interactor2; # /; 
    my @interactor2 = ();
    foreach my $node (@tmp) {
	foreach (split /;/, $node) { #/) {
	    push @interactor2, $_;
	}
    }
    $all_connections++;

    foreach my $interactor1 (@interactor1) {
	foreach my $interactor2 (@interactor2) {

	    # store interaction in both directions
	    $connections{$interactor1}{$interactor2}++; 
	    $connections{$interactor2}{$interactor1}++; 

	    # if both genes are amongst the DEGs, keep them in another hash:
	    if (defined $genes{$interactor1}
		and
		defined $genes{$interactor2}) {
		$deg_connections++;
		$interactions{$interactor1}{$interactor2}++;
		$interactions{$interactor2}{$interactor1}++;
	    }	    
	} 
    }
}
close IN;

my $all_genes = scalar keys %interactions;
my $network_genes = scalar keys %connections;
print STDERR "Found $all_genes genes with $deg_connections connections amongst $network_genes nodes and $all_connections edges from network\n";

if ($network_genes > $all_genes) {
    print STDERR "Nodes not in DEG list:\n";

    foreach my $node (sort keys %connections) {
	unless (defined $interactions{$node}) {
	    print STDERR "$node\n";
	}
    }
}

# make sure all hub/bottleneck genes were found:
foreach my $hub (@hubs) {
    unless (defined $interactions{$hub}) {
	die "Gene '$hub' not found in network!\n";
    }
}

# Step 3: get numbers of directly connected genes:

print "Extension\thubs\tDEG genes\tpercentage of all\tpercentage of affected\tGenes\n";

# create a hash that contains the connected genes as keys and the hub/bottleneck genes as values
# this way we can capture genes linked to multiple hub/bottleneck genes
my %lists = ();
foreach my $hub (@hubs) {
    my @link1 = (sort keys %{$interactions{$hub}});
    foreach (@link1) {
	$lists{$_}{$hub}++;
    }
}

# now turn this list around by using the hub/bottleneck genes as keys
# and the genes linked to them as values:
my %list2 = ();
foreach my $gene (sort keys %lists) {
    my $hubs = join ", ", sort keys %{$lists{$gene}};
    $list2{$hubs}{$gene}++;
}

# finally, for each sub-list print out the details
# (and capture information for next link level):
my $links = 0;
my $sum = 0;

# get number of affected genes first
# so we can calcalate the percentage in the next run
foreach my $hubs (sort @hubs) {
    my $num = 1;
    $sum += $num;
}

my $perc_sum = 0;
foreach my $hubs (sort @hubs) {
    my $num = 1;
    my $perc = sprintf "%.1f", $num / $all_genes * 100;
    my $genes = $hubs;
    my $perc2 = sprintf "%.1f", $num / $sum * 100;
    $perc_sum += $num / $sum * 100;
    print "$links\t$hubs\t$num\t$perc\t$perc2\t$genes\n";
}
my $perc = sprintf "%.1f", $sum / $all_genes * 100;
my $total_genes = join ', ', @hubs;
print "$links\ttotal\t$sum\t$perc\t$perc_sum\t$total_genes\n";

# determine the number of genes that are not linked to hub/bottleneck genes:
my %missing = ();
foreach (sort keys %interactions) {
    unless (defined $hubs{$_}) {
	$missing{$_}++;
    }
}
my $missing = scalar keys %missing;
my $missing_genes = join ', ', sort keys %missing;
$perc = sprintf "%.1f", $missing / $all_genes * 100;
print "$links\tunconnected\t$missing\t$perc\t-\t$missing_genes\n";

$links++;
print STDERR "Extend to $links connections ($sum genes)\n";

$sum = 0;
# get number of affected genes first
# so we can calcalate the percentage in the next run
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    $sum += $num;
}

my %list_next = ();
$perc_sum = 0;
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    foreach (keys %{$list2{$hubs}}) {
	$list_next{$_}++;
    }
    my $perc = sprintf "%.1f", $num / $all_genes * 100;
    my $perc2 = sprintf "%.1f", $num / $sum * 100;
    my $genes = join ', ', sort keys %{$list2{$hubs}};
    $perc_sum += $num / $sum * 100;
    print "$links\t$hubs\t$num\t$perc\t$perc2\t$genes\n";
}
$perc = sprintf "%.1f", $sum / $all_genes * 100;
$total_genes = join ', ', sort keys %list_next;
print "$links\ttotal\t$sum\t$perc\t$perc_sum\t$total_genes\n";

# %list_next now contains all the differentially expressed genes linked at this level to the hub/bottleneck genes

# determine the number of genes that are not linked to hub/bottleneck genes:
%missing = ();
foreach (sort keys %interactions) {
    unless (defined $list_next{$_}) {
	$missing{$_}++;
    }
}
$missing = scalar keys %missing;
$missing_genes = join ', ', sort keys %missing;
$perc = sprintf "%.1f", $missing / $all_genes * 100;
print "$links\tunconnected\t$missing\t$perc\t-\t$missing_genes\n";


# Step 4: get number of genes connected via up to 2 links:

# (something more elegant could have been done via recursion,
# but for a limited number of steps this explicit approach was easier)

my @list_next = sort keys %list_next;
my $next = scalar @list_next;
$links++;
print STDERR "Extend to $links connections ($next genes)\n";

%lists = ();
foreach my $hub (@hubs) {
    # 1st degree neighbours
    my @link1 = (sort keys %{$interactions{$hub}});
    foreach my $hub2 ($hub, @link1) {
	# 2nd degree neighbours
	$lists{$hub2}{$hub}++;
	foreach (sort keys %{$interactions{$hub2}}) {
	    $lists{$_}{$hub}++;
	}
    }
}

%list2 = ();
foreach my $gene (sort keys %lists) {
    my $hubs = join ", ", sort keys %{$lists{$gene}};
    $list2{$hubs}{$gene}++;
}

# convert list for output sorted by hub combinations:

# get number of affected genes first
# so we can calcalate the percentage in the next run
$sum = 0;
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    $sum += $num;
}

%list_next = ();
$perc_sum = 0;
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    foreach (keys %{$list2{$hubs}}) {
	$list_next{$_}++;
    }
    my $perc = sprintf "%.1f", $num / $all_genes * 100;
    my $perc2 = sprintf "%.1f", $num / $sum * 100;
    my $genes = join ', ', sort keys %{$list2{$hubs}};
    $perc_sum += $num / $sum * 100;
    print "$links\t$hubs\t$num\t$perc\t$perc2\t$genes\n";
}
$perc = sprintf "%.1f", $sum / $all_genes * 100;
$total_genes = join ', ', sort keys %list_next;
print "$links\ttotal\t$sum\t$perc\t$perc_sum\t$total_genes\n";

%missing = ();
foreach (sort keys %interactions) {
    unless (defined $list_next{$_}) {
        $missing{$_}++;
    }
}
$missing = scalar keys %missing;
$missing_genes = join ', ', sort keys %missing;
$perc = sprintf "%.1f", $missing / $all_genes * 100;
print "$links\tunconnected\t$missing\t$perc\t-\t$missing_genes\n";


# Step 5: get number of genes connected via up to 3 links:

@list_next = sort keys %list_next;
$next = scalar @list_next;
$links++;
print STDERR "Extend to $links connections ($next genes)\n";
%lists = ();
foreach my $hub (@hubs) {
    # 1st degree neighbours 
    my @link1 = (sort keys %{$interactions{$hub}});
    my %lst = ();
    foreach my $hub2 ($hub, @link1) {
	# 2nd degree neighbours
	$lst{$hub2}++;
        foreach (sort keys %{$interactions{$hub2}}) {
            $lst{$_}++;
        }
    }
    foreach my $hub3 (sort keys %lst) {
	# 3rd degree neighbours
	$lists{$hub3}{$hub}++;
	foreach (sort keys %{$interactions{$hub3}}) {
	    $lists{$_}{$hub}++;
	}
    }
}

%list2 = ();
foreach my $gene (sort keys %lists) {
    my $hubs = join ", ", sort keys %{$lists{$gene}};
    $list2{$hubs}{$gene}++;
}


$sum = 0;
# get number of affected genes first
# so we can calcalate the percentage in the next run
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
   $sum += $num;
}

%list_next = ();
$perc_sum = 0;
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    foreach (keys %{$list2{$hubs}}) {
        $list_next{$_}++;
    }
    my $perc = sprintf "%.1f", $num / $all_genes * 100;
    my $perc2 = sprintf "%.1f", $num / $sum * 100;
    my $genes = join ', ', sort keys %{$list2{$hubs}};
    $perc_sum += $num / $sum * 100;
    print "$links\t$hubs\t$num\t$perc\t$perc2\t$genes\n";
}
$perc = sprintf "%.1f", $sum / $all_genes * 100;
$total_genes = join ', ', sort keys %list_next;
print "$links\ttotal\t$sum\t$perc\t$perc_sum\t$total_genes\n";

%missing = ();
foreach (sort keys %interactions) {
    unless (defined $list_next{$_}) {
        $missing{$_}++;
    }
}
$missing = scalar keys %missing;
$missing_genes = join ', ', sort keys %missing;
$perc = sprintf "%.1f", $missing / $all_genes * 100;
print "$links\tunconnected\t$missing\t$perc\t-\t$missing_genes\n";


# Step 6: get number of genes connected via up to 4 links:

@list_next = sort keys %list_next;
$next = scalar @list_next;
$links++;
print STDERR "Extend to $links connections ($next genes)\n";
%lists = ();
foreach my $hub (@hubs) {
    # 1st degree neighbours
    my @link1 = (sort keys %{$interactions{$hub}});
    my %lst2 = ();
    foreach my $hub2 ($hub, @link1) {	
	# 2nd degree neighbours
	$lst2{$hub2}++;
        foreach (sort keys %{$interactions{$hub2}}) {
            $lst2{$_}++;
        }
    }

    my %lst = ();
    foreach my $hub2 (sort keys %lst2) {
	# 3rd degree neighbours
        $lst{$hub2}++;
        foreach (sort keys %{$interactions{$hub2}}) {
            $lst{$_}++;
        }
    }

    foreach my $hub3 (sort keys %lst) {
	# 4th degree neighbours
	$lists{$hub3}{$hub}++;
	foreach (sort keys %{$interactions{$hub3}}) {
	    $lists{$_}{$hub}++;
	}
    }
}

%list2 = ();
foreach my $gene (sort keys %lists) {
    my $hubs = join ", ", sort keys %{$lists{$gene}};
    $list2{$hubs}{$gene}++;
}

$sum = 0;
# get number of affected genes first
# so we can calcalate the percentage in the next run
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    $sum += $num;
}

%list_next = ();
$perc_sum = 0;
foreach my $hubs (sort keys %list2) {
    my $num = scalar keys %{$list2{$hubs}};
    foreach (keys %{$list2{$hubs}}) {
        $list_next{$_}++;
    }
    my $perc = sprintf "%.1f", $num / $all_genes * 100;
    my $perc2 = sprintf "%.1f", $num / $sum * 100;
    my $genes = join ', ', sort keys %{$list2{$hubs}};
    $perc_sum += $num / $sum * 100;
    print "$links\t$hubs\t$num\t$perc\t$perc2\t$genes\n";
}
$perc = sprintf "%.1f", $sum / $all_genes * 100;
$total_genes = join ', ', sort keys %list_next;
print "$links\ttotal\t$sum\t$perc\t$perc_sum\t$total_genes\n";
%missing = ();
foreach (sort keys %interactions) {
    unless (defined $list_next{$_}) {
        $missing{$_}++;
    }
}
$missing = scalar keys %missing;
$missing_genes = join ', ', sort keys %missing;
$perc = sprintf "%.1f", $missing / $all_genes * 100;
print "$links\tunconnected\t$missing\t$perc\t-\t$missing_genes\n";
