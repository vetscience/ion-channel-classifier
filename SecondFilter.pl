#SecondFilter.pl @ Bahiyah Nor, Aug 18th 2015
#
# Usage: SecondFilter.pl -i <BLAST.tabdelimited.output>
# Other option: -h For Help information
# 
# Functions similarly to FirstFilter.pl, excepts that it
# processes results from BLASTp against MuSICC Ion Channel
# and Aquaporin Database.
# Produces a report containing the Sequence Ids and the top hit
# Will be useful for sequence grouping process.
# Output: SecondFilter.Report.txt
use strict;
use warnings;

my $blastfile = "";
my $matchArg = 0;

foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i')
  {
    if($i+1 < @ARGV)
    {
      $blastfile = $ARGV[$i+1];
      open TMP, $blastfile or die "Error encountered. $blastfile cannot be loaded.\n";
      close TMP;
      $matchArg = 1;
    }
    else
    {
      print "\n***!!! Error encountered. Invalid filename input.\n\n";
    }
  }
  elsif($ARGV[$i] eq '-h')
  {
    print "#SecondFilter.pl @ Bahiyah Nor, Aug 18th 2015\n\n";
    print "# Usage: SecondFilter.pl -i <BLAST.tabdelimited.output>\n";
    print "# Other option: -h For Help information\n\n";
    print "# Functions similarly to FirstFilter.pl, excepts that it\n";
    print "# processes results from BLASTp against MuSICC Ion Channel\n";
    print "# and Aquaporin Database.\n";
    print "# Produces a report containing the Sequence Ids and the top hit\n";
    print "# Will be useful for sequence grouping process.\n";
    die "# Output: SecondFilter.Report.txt\n\n";
  }
}

if(@ARGV == 0 or $matchArg == 0)
{
  print "#SecondFilter.pl @ Bahiyah Nor, Aug 18th 2015\n\n";
  print "# Usage: SecondFilter.pl -i <BLAST.tabdelimited.output>\n";
  print "# Other option: -h For Help information\n\n";
  print "# Functions similarly to FirstFilter.pl, excepts that it\n";
  print "# processes results from BLASTp against MuSICC Ion Channel\n";
  print "# and Aquaporin Database.\n";
  print "# Produces a report containing the Sequence Ids and the top hit\n";
  print "# Will be useful for sequence grouping process.\n";
  die "# Output: SecondFilter.Report.txt\n\n";
}

#Load the KEGG Ion Channel and Aquaporin K-terms
my $ktermfile = "Database/KEGG.Kterms.txt";
open KEGG, $ktermfile or die "Cannot load $ktermfile\n";
my %KTERMS;
my $kegg;
my $preFam;

while($kegg=<KEGG>)
{
  chomp($kegg);
  my @parts = split(/\s+/,$kegg);
  if($parts[0] eq 'B')
  {
    $preFam = $parts[1];
  }
  else
  {
    $KTERMS{$parts[1]} = $preFam;
  }
}

close KEGG;

#Load the blast result file, pick the best hit and write the report
my $output = "Output/SecondFilter.Report.txt";
if(-e $output)
{
  unlink $output;
}
open WRITER, ">>$output" or die "Cannot create $output.\n";
my $writer;

open BLAST, $blastfile or die "Cannot load $blastfile.\n";
my $blast;
my %RESULTS;

while($blast=<BLAST>)
{
  chomp($blast);
  my @cells = split(/\s+/,$blast);
  
  #Check the file format (BLASTp tab-delimited output has 12 columns)
  #Program will terminate if this criterion is not met
  if(@cells != 12)
  {
    die "The format output of $blastfile is not recognized by the program.\nMuSICC requires tab-delimited output from BLASTp.\n";
  }
  
  if(!exists $RESULTS{$cells[0]})
  {
    my $ko = substr($cells[1],0,6);
    $RESULTS{$cells[0]} = $KTERMS{$ko};
  }
}

close BLAST;

foreach my $seq (keys %RESULTS)
{
  $writer = $seq . "\t" . $RESULTS{$seq} . "\n";
  print WRITER $writer;
}

close WRITER;
