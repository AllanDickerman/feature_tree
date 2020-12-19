use Sequence_Alignment;
use Data::Dumper qw(Dumper);
 
print Dumper \@ARGV;

if (scalar(@ARGV) < 2) {
    print "Usage: trim_ends_and_delete_bad_seqs.pl alignment_file end_dendisty_threshold [row_density_thresh]\n";
    exit(1)
}

my $file = $ARGV[0];
my $trim_threshold = $ARGV[1];
my $row_threshold = $trim_threshold;
$row_thrshold = $ARGV[2] if @ARGV > 2;
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

$al->end_trim($trim_threshold);

$al->delete_gappy_seqs($row_threshold);

open $fh, ">", $file."_trimmed.afa";
$al->write_fasta($fh)
