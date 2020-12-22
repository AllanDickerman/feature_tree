use Sequence_Alignment;
use Data::Dumper qw(Dumper);
 
print Dumper \@ARGV;

if (scalar(@ARGV) < 2) {
    exit(1)
}

my $file = $ARGV[0];
my $trim_threshold = $ARGV[1];
my $fh;
unless (open($fh, $file)) {
    die "Problem, could not open $file\n"
}
$al = new Sequence_Alignment($fh);
print "al = $al\n";
print "Number of taxa = ", $al->get_ntaxa(), "\n";
print "Length of alignment = ", $al->get_length(), "\n";
print "Is aligned = ", $al->is_aligned(), "\n";
print "al = $al\n";

$al->end_trim(0.5);
$al->write_fasta($file."_end_trimmed.afa")
