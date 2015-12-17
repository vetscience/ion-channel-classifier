#PseAAC.DataPreProcess.pl @ Bahiyah Nor, Aug 18th 2015
#
# Usage: PseAAC.DataPreProcess.pl -i <fastafile>
# Other option: -o Output.fa  The output filename
# The program converts the fasta file into format recognized by
# BuildPseAAC.pl.
# This is necessary to ensure that the PseAAC for classification
# is computed correctly.
# Ensure that the sequence is in one line, instead of broken in
# different lines.
# Output: FormattedFasta.fa
use strict;
use warnings;

my $input = "";
my $inputFlag = 0;
my $output = "Output/FormattedFasta.fa";
my $matchArg = 0;

foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i')
  {
    if($i+1 < @ARGV)
    {
      $input = $ARGV[$i+1];
      open TMP, $input or die "Cannot load $input. Please ensure that the file exists.\n";
      close TMP;
      $matchArg = 1;
      $inputFlag = 1;
    }
    else
    {
      print "**!! Error Encountered. Invalid file provided.\n";
    }
  }
  elsif($ARGV[$i] eq '-o')
  {
    if($i+1 < @ARGV)
    {
      $output = $ARGV[$i+1];
      $matchArg = 1;
    }
    else
    {
      print "No output filename is provided. The default output file will be used.\n";
      $matchArg = 1;
    }
  }
}

if(@ARGV == 0 or $inputFlag == 0 or $matchArg == 0)
{
  print "#PseAAC.DataPreProcess.pl @ Bahiyah Nor, Aug 18th 2015\n#\n";
  print "# Usage: PseAAC.DataPreProcess.pl -i <fastafile>\n";
  print "# Other option: -o Output.fa  The output filename\n";
  print "# The program converts the fasta file into format recognized by\n";
  print "# BuildPseAAC.pl.\n";
  print "# This is necessary to ensure that the PseAAC for classification\n";
  print "# is computed correctly.\n";
  print "# Ensure that the sequence is in one line, instead of broken in\n";
  print "# different lines.\n";
  die "# Default Output File: FormattedFasta.fa\n";
}

#Load the file and reformat the file
open SEQ, $input or die "Cannot load $input.\n";
my $seq;
my $read = "";
my $id = "";

if(-e $output)
{
  unlink $output;
}
open WRITER, ">>$output" or die "Cannot create $output.\n";
my $writer;

while($seq=<SEQ>)
{
  chomp($seq);
  if(substr($seq,0,1) eq '>') #The fasta header
  {
    if(length($id) > 0)
    {
      $writer = $id . "\n" . $read . "\n";
      print WRITER $writer;
      $read = "";
    }
    $id = $seq;
    $read = "";
  }
  else
  {
    if(length($seq) > 0)
    {
      $read = $read . $seq;
    }
  }

}

if(length($read) > 0)
{
  $writer = $id . "\n" . $read;
  print WRITER $writer;
}

close SEQ;
close WRITER;
