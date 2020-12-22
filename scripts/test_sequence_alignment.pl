use Sequence_Alignment;
use Data::Dumper qw(Dumper);
 
print Dumper \@ARGV;

my $file = $ARGV[0];

print STDERR "In test: file = $file\n";
my $fh;
unless (open($fh, $file)) {
    die "Problem, could not open $file\n"
}
print STDERR "In test: fh = $fh\n";
$al = new Sequence_Alignment($fh);
print "al = $al\n";
print "Number of taxa = ", $al->get_ntaxa(), "\n";
print "Length of alignment = ", $al->get_length(), "\n";
print "Is aligned = ", $al->is_aligned(), "\n";
print "al = $al\n";
my $fh;
open $fh, ">", "test.phy";
print "fh = $fh\n";
$al->write_phylip($fh)
