#FirstFilter.pl @ Bahiyah Nor, Aug 17th 2015
#
# Usage: FirstFilter.pl -i <blastp.tabdelimited.results>
# [Optional parameters: 1)-e <E-value threshold> 2)-h for help]
#
# ** Step 2 of MuSICC: Initial screening for homologous sequences
# Screen the results from BLASTp results to the KEGG database
# Select sequences homologous to ion channels based on K-terms
# limited by the e-value threshold
# Default e-value threshold: 1E-15
# Two files generated: 1. Text file containing headers for sequence retrieval
#                      2. Report file with sequence ids, hits and e-values
use strict;
use warnings;

my $blastres = "";
my $evalue = 1E-15;
my $inputFlag = 0;
my $matchARG = 0;

#Grab all the input arguments
foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i') #The blast results file is passed in the next argument
  {
    if($i+1 < @ARGV)
    {
      $blastres = $ARGV[$i+1];
      open TMP, $blastres or die "Cannot load $blastres, please make sure that it is the right file.\n";
      close TMP;
      $inputFlag = 1;
      $matchARG = 1;
    }
    else
    {
      print "Error encountered. The BLAST result file was not specified.\n";
    }
  }
  if($ARGV[$i] eq '-e') #User-defined e-value threshold
  {
    if($i+1 < @ARGV)
    {
      $evalue = $ARGV[$i+1];
      $matchARG = 1;
    }
    else
    {
      print "The '-e' option was chosen without a specified threshold. The default value will be used.\n";
    }
  }
  if($ARGV[$i] eq '-h') #Help requested
  {
    print "\nFirstFilter.pl @ MuSICC First Screening Process\n\nUsage: FirstFilter.pl -i <blastp.tabdelimited.results>\n";
    print "---------------------------------------------------------------\n";
    print "\tOther options:\t-e <E-value threshold>\n\t\t\t-h For help\n";
    print "\nFirstFilter.pl screens for homologous sequences base on the blastp results.\n";
    die "FirstFilter.pl is terminated.\n\n";
  }
}

#If no arguments were passed
if(@ARGV == 0 or $inputFlag == 0 or $matchARG == 0)
{
  print "\nFirstFilter.pl @ MuSICC First Screening Process\n\nUsage: FirstFilter.pl -i <blastp.tabdelimited.results>\n";
  print "---------------------------------------------------------------\n";
  print "\tOther options:\t-e <E-value threshold>\n\t\t\t-h For help\n";
  print "\nFirstFilter.pl screens for homologous sequences base on the blastp results.\n";
  die "FirstFilter.pl is terminated.\n\n";
}

#Load the KEGG Ion Channel and  Gene List
#If the gene list from KEGG has been updated, rename the new list with the same file name
#or make the changes here
my $GeneFile = "Database/KEGG.GeneList.txt";
open GENELIST, $GeneFile or die "Cannot load $GeneFile, please ensure that the file exists.\n";
my %GeneList;
my $gene;

print "\tLoading the KEGG Gene List from $GeneFile...\n";
while($gene=<GENELIST>)
{
  chomp($gene);
  my ($kterm,$id) = split(/\s+/,$gene);
  $GeneList{$id} = $kterm;
}

close GENELIST;

#Load the BLASTp results
#The program will check if the results is the right format (tab-delimited)
open BLAST, $blastres or die "Cannot load $blastres. Please ensure that the file exists.\n";
my %Headers;
my $blast;

print "\tProcessing $blastres...\n";
while($blast=<BLAST>)
{
  chomp($blast);
  my @cells = split(/\s+/,$blast);
  
  #Check the file format (BLASTp tab-delimited output has 12 columns)
  #Program will terminate if this criterion is not met
  if(@cells != 12)
  {
    die "The format output of $blastres is not recognized by the program.\nMuSICC requires tab-delimited output from BLASTp.\n";
  }

  #Check if the sequence hit an ion channel gene below the e-value threshold
  if(!exists $Headers{$cells[0]}) #The sequence has not been added to the Headers list
  {
    if($cells[10] < $evalue)
    {
      if(exists $GeneList{$cells[1]}) #If the gene hit an ion channel
      {
        $Headers{$cells[0]}{'Evalue'} = $cells[10];
        $Headers{$cells[0]}{'Hit'} = $cells[1];
      }
    }
  }  
}

close BLAST;

#Print output to files
#Two files will be generated: 1. Text file containing the headers only for sequence retrieval (FirstFilter.Headers.txt)
#                             2. Report file containing Sequence IDs, Hit and Evalue (FirstFilter.Report.txt)
#Existing files of the same name in this directory will be removed.
my $report = "Output/FirstFilter.Report.txt";
my $head = "Output/FirstFilter.Headers.txt";
if(-e $report)
{
  unlink $report;
}
if(-e $head)
{
  unlink $head;
}
open WRITER_REP, ">>$report" or die "Error encountered. Cannot create $report.\n";
open WRITER_HEAD, ">>$head" or die "Error encountered. Cannot create $head.\n";

print "\tPrinting results to $report and $head...\n";

foreach my $seqId (keys %Headers)
{
  my $head_line = $seqId . "\n";
  print WRITER_HEAD $head_line;

  my $rep_line = $seqId . "\t" . $Headers{$seqId}{'Hit'} . "\t" . $Headers{$seqId}{'Evalue'} . "\n";
  print WRITER_REP $rep_line;
}

close WRITER_REP;
close WRITER_HEAD;

print "\tEnd of filtering. Results can be found in $report and $head\n";
