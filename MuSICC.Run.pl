#!/usr/bin/perl
#MuSICC.Run.pl @ Bahiyah Nor, Aug 17th 2015
#This scripts runs the whole MuSICC pipeline using default threshold values
#To run parts of the pipeline, the perl scripts are available in the Scripts folder
#
# Usage: MuSICC.Run.sh -i <protein fasta file>
# Options: -n  The number of threads to be used. Default = 1 for all the processes.
#              The threads will be divided between processes
#
#          -b  BLASTp tab-delimited output file. 
#              If this file is provided, the first BLASTp process will be skipped.
#         
#		   -k  The path to the KEGG database if -b option not provided 
#
#          -t  The path to the training set blast database
#
#          -h  Help 
#
# The initial Protein FASTA file is required by MuSICC.
use strict;
use warnings;
use POSIX;

my $fasta = "";
my $num_thread = 1;
my $blast = "";
my $ionDB = "";
my $keggdatabase = "";
#Pass the command line arguments
foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i')
  {
    $fasta = $ARGV[$i+1];
    open TMP, $fasta or die "Cannot load $fasta. Please ensure that the file exists.\n";
    close TMP;
  }
  if($ARGV[$i] eq '-n')
  {
    $num_thread = $ARGV[$i+1];
  }
  if($ARGV[$i] eq '-b')
  {
    $blast = $ARGV[$i+1];
  }
  if($ARGV[$i] eq '-k')
  {
    $keggdatabase = $ARGV[$i+1];
  }
  if($ARGV[$i] eq '-t')
  {
    $ionDB = $ARGV[$i+1];
  }
  if($ARGV[$i] eq '-h')
  {
    system("cat README.txt");
    exit();
  }
}

if(@ARGV == 0)
{
  message();
}

my $launchTime = strftime "%Y%m%d_%H%M%S", localtime;
my $logFile = "logs/" . $launchTime . ".log";
my $errorfile = "logs/" . $launchTime . ".error.log";
open SYSLOG, ">>$logFile" or die "Cannot create $logFile\n";
my $syslogging = "MuSICC run for :" . $fasta . " at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
print "\n" . $syslogging;
print "------------------------------------------------------\n";
print SYSLOG $syslogging;

#Step 1-3: BLASTp against the KEGG DATABASE
# KEGG database location now on command line 15/12/2015, Ross
# my $keggdatabase = "/home/Public/Ps1/BlastDb/Kegg2015Feb/bah_ion"; # modify this to the location of the kegg database

if($blast eq "") #BLASTp has not been run, results not provided
{
  $syslogging = "Step 1: Start of BLASTp on $fasta " . localtime() . "\n";
  print $syslogging;
  print SYSLOG $syslogging;
  system("blastp -query $fasta -db $keggdatabase -num_threads $num_thread -evalue 1e-15 -outfmt 6 -out Output/FirstBLASTp.csv 2> $errorfile");
  if($? != 0)
  {
    $syslogging = "\tERROR: BLASTp against KEGG database, error at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
    print SYSLOG $syslogging;
    die $syslogging;
  }
  exit();
  $syslogging = "\tEnd of BLASTp " . strftime("%Y %b %a %H:%M:%S", localtime) . ".\nStep 2: Start of first filtering process\n";
  print $syslogging;
  print SYSLOG $syslogging;
  system("perl FirstFilter.pl -i Output/FirstBLASTp.csv 2> $errorfile");
  if($? != 0)
  {
    $syslogging = "\tERROR: FirstFilter.pl encounters an error at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
    print SYSLOG $syslogging;
    die $syslogging;
  }
  $syslogging = "Step 3: Retrieving filtered sequences at ." . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
  print $syslogging;
  print SYSLOG $syslogging;
  system("perl RetrieveSeqs.pl -f $fasta Output/FirstFilter.Headers.txt > Output/FirstFilter.Sequences.fa 2> $errorfile");
  if($? != 0)
  {
    $syslogging = "\tERROR: RetrieveSeqs.pl encounters an error at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
    print SYSLOG $syslogging;
    die $syslogging;
  }
}
else
{
  #This is Exception 1: BLASTp results provided
  $syslogging = "BLASTp result file provided.\nSkipping Steps 1 and 2...\n";
  print $syslogging;
  print SYSLOG $syslogging;
  $syslogging = "Extracting Sequence IDs from $blast " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
  print $syslogging;
  print SYSLOG $syslogging;
  my $out = "Output/Exception1.Headers.txt";
  system("awk '{print \$1}' $blast | uniq > Output/Exception1.Headers.txt 2> $errorfile");
  if($? != 0)
  {
    $syslogging = "\tERROR: MuSICC encounters error while extracting sequence IDs at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
    print SYSLOG $syslogging;
    die $syslogging;
  }
  $syslogging = "Step 3: Retrieving filtered sequences at ." . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
  print $syslogging;
  print SYSLOG $syslogging;
  system("perl RetrieveSeqs.pl -f $fasta Output/Exception1.Headers.txt > Output/FirstFilter.Sequences.fa");
  if($? != 0)
  {
    $syslogging = "\tERROR: RetrieveSeqs.pl encounters an error at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
    print SYSLOG $syslogging;
    die $syslogging;
  }
}

#Step 4-5: BLASTp to MuSICC Ion Channel and Aquaporin Database
# Training set database location now on command line 15/12/2015, Ross 
# my $ionDB = "Database/MuSICC.IonChannels.fa";
$syslogging = "Step 4: BLASTp against MuSICC Ion Channel and Aquaporin Database ". strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("blastp -query Output/FirstFilter.Sequences.fa -db $ionDB -evalue 1e-45 -outfmt 6 -out Output/SecondBLASTp.csv -num_threads $num_thread 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: BLASTp against MuSICC Ion Channel and Aquaporin Database encounters an error at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}
system("perl SecondFilter.pl -i Output/SecondBLASTp.csv 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while filtering BLASTp results at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}
$syslogging = "Step 5: Retrieving filtered sequences " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("awk '{print \$1}' Output/SecondBLASTp.csv | uniq > Output/SecondBLASTp.Headers.txt 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while retrieving sequences at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}
system("perl RetrieveSeqs.pl -f Output/FirstFilter.Sequences.fa Output/SecondBLASTp.Headers.txt > Output/FilteredSequences.fa 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while filtering sequences at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

#Step 6: Conserved Domain Prediction
$syslogging = "Step 6: Predicting the conserved domains " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("interproscan.h -i Output/FilteredSequences.fa -f tsv --goterms -dp -o Output/ipro.tsv 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while predicting conserved domains at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

#Step 7: Transmembrane Domain Prediction using TMHMM
$syslogging = "Step 7: Predicting the transmembrane domain prediction " . strftime("%Y %b %a %H:%M:%S",localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("tmhmm -short Output/FilteredSequences.fa > Output/TM.Predicted.txt 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while predicting transmembrane domains at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

#Step 8: Grouping of the sequences base on the filtration criteria
$syslogging = "Step 8: Grouping the sequences based on the filtration criteria " . strftime("%Y %b %a %H:%M:%S",localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("perl SequenceGrouping.pl -i Output/SecondFilter.Report.txt -c Output/ipro.tsv -t Output/TM.Predicted.txt 2> $errorfile"); 
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while forming the groups at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

#Step 9: Compute PseAAC
$syslogging = "Step 9: Computing the pseudo amino acid composition " . strftime("%Y %b %a %H:%M:%S",localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("perl PseAAC.DataPreProcess.pl -i Output/FilteredSequences.fa 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while formatting the sequences at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}
system("perl BuildPseAAC.pl -i Scripts/Output/FormattedFasta.fa 2> $errorfile");
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while computing the PseAAC at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

#Step 10: Classification via SVM
$syslogging = "Step 10: Classifying the ion channels " . strftime("%Y %b %a %H:%M:%S",localtime) . "\n";
print $syslogging;
print SYSLOG $syslogging;
system("Rscript SVMClassifier.r 2> $errorfile"); #Results can be found in MuSICC.Results.date.txt
if($? != 0)
{
  $syslogging = "\tERROR: MuSICC encounters an error while classifying ion channels at " . strftime("%Y %b %a %H:%M:%S", localtime) . "\n\tError message can be found in $errorfile.\n\n";
  print SYSLOG $syslogging;
  die $syslogging;
}

sub message {

  print "\nMuSICC.Run.pl @ Multi-Screening Ion Channel Classifier\n";
  print "This program runs the whole MuSICC pipeline.\n";
  print "\nUsage: MuSICC.Run.sh -i <protein fasta file>\n";
  print "Options:\t-n\tThe number of threads to be used. Default = 1 for all the processes.\n";
  print "\t\t\tThe threads will be divided between processes\n";
  print "\t\t-b\tBLASTp tab-delimited output file.\n";
  print "\t\t\tIf this file is provided, the first BLASTp process will be skipped.\n";
  print "\t\t\tThe initial protein FASTA file has to be provided\n";
  print "\t\t-k\tThe path to the KEGG database if -b option not provided.\n";
  print "\t\t-t\tThe path to the training set blast database.\n";    
  print "\t\t-h\tHelp\n";
  print "\n------------------------------------------------------\n";
  print "\n MuSICC is a multi-screening ion channel classifier.\n";
  print " The program screens for potential ion channel and aquaporin sequences in the FASTA file provided\n";
  print " based on three defined criteria:\n";
  print "\t1. Homology to known ion channels/aquaporins\n";
  print "\t2. Presence of known ion channels/aquaporins conserved domains\n";
  print "\t3. Transmembrane domain counts within the expected ranges of known ion channels/aquaporins\n";
  print " MuSICC uses SVM-based classifier to classify the unknown sequences\n";
  print " There are 48 ion channel subfamilies and 1 aquaporin subfamily included in the classifier.\n"; 
  print "\n** MuSICC requires the following Bioinformatics tools installed for a complete run:\n";
  print "\t1. BLASTp\n\t2. KEGG database\n\t3. InterProScan version 5.7.48\n\t4. TMHMM 2.0\n";
  print " Please refer to the README.txt file for other options in cases where the tools are not available.\n";
  die "------------------------------------------------------\n";

}
