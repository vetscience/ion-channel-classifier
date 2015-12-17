#CompileResults.pl @ Bahiyah Nor, Sept 15th 2015
#
# Usage: CompileResults.pl -s <SVMClassification.res> -g <sequencegrouping>
# Options: -o <Output.tag> All the outputs from this program will be tagged as such
#	   -h For help
# Compile the results of MuSICC run
# Organise the classification made through the SVM classifier
# with the groupings made
# Present the results in a table with sequence id, classification, probability value
# and the grouping (refer to manual regarding grouping)
# Outputs: 1. MuSICC.[launchtime].Classification.txt 2. MuSICC.[launchtime].ChannelCount.txt 
use strict;
use warnings;
use POSIX;

my $svmresfile;
my $svmFlag = 0;
my $seqgroupfile;
my $groupFlag = 0;
my $launchTime = strftime "%Y%m%d_%H%M%S", localtime;
my $output1 = "MuSICC." . $launchTime . ".Classification.txt";
my $output2 = "MuSICC." . $launchTime . ".ChannelCount.txt";

foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-s')
  {
    $svmresfile = $ARGV[$i+1] or message();
    $svmFlag = 1;
  }
  if($ARGV[$i] eq '-g')
  {
    $seqgroupfile = $ARGV[$i+1] or message();
    $groupFlag = 1;
  }
  if($ARGV[$i] eq '-o')
  {
    if($i + 1 > $#ARGV)
    {
      print "Output tag was inserted but none given. Default filenames will be used.\n";
    }
    else
    {
      $output1 = $ARGV[$i+1] . ".Classification.txt";
      $output2 = $ARGV[$i+1] . ".ChannelCount.txt";
    }
  }
  if($ARGV[$i] eq '-h')
  {
    message();
  }
}

if($svmFlag == 0 or $groupFlag == 0 or @ARGV == 0) 
{
  message();
}

#Load the Ion Channel subfamily to count
my @ionFam = ("GABA","GLYR","ANGL","ACHN","SERO","ZNAC","GLUT_NNMDA","GLUT_NMDA","EPIT","ASIC","ATPG","RYAN","IP3R","SCNA","SCNB","CACNA1","CACNA2","CACNB","CACNG","SHAKER","SHAB","SHAW","SHAL","ISK","KCNF","KCNG","KCNH","KCNK","KCNQ","KCNS","KCNV","KCNM","KCNN","KCNT","IRCHN","CNG","HCN","TRPC","TRPV","TRPM","TRPA","TRPP","TRPML","CAT2P","CLCN","NSCCN","CICH","CACC","AQUA");
my %ionCount;
foreach my $chan (0 .. $#ionFam)
{
  $ionCount{$ionFam[$chan]} = 0;
}

#load the results
open SVM, $svmresfile or die "Cannot load $svmresfile. Please ensure that this is the right file.\n";
my $svm;
my %SVMclass;
my @SequenceID;
my $c = 0;
my $totalChan = 0;

while($svm = <SVM>)
{
  chomp($svm);
  if($c > 0)
  {
    my @cells = split(/\s+/, $svm);
    my $seqId = substr($cells[1],2,length($cells[1])-3);
    push(@SequenceID, $seqId);
    my $classified = substr($cells[4],1,length($cells[4])-2);
    $SVMclass{$seqId}{'Class'} = $classified;
    $SVMclass{$seqId}{'pval'} = $cells[3];
    $ionCount{$classified} = $ionCount{$classified} + 1;
    $totalChan = $totalChan + 1; 
  }
  else
  {
    my @tmp = split(/\s+/,$svm);
    if(@tmp != 4)
    {
      die "The file provided is not the result from MuSICC SVM classifier. Please ensure that you have the right file.\n";
    }
  }
  $c = $c + 1;
}

close SVM;

open GROUP, $seqgroupfile or die "Cannot load $seqgroupfile. Please ensure that this is the right file.\n";
my $group;
my $line = 0;

while($group=<GROUP>)
{
  chomp($group);
  if($line > 0)
  {
    my @gcells = split(/\s+/,$group);
    if(@gcells != 5)
    {
      die "$seqgroupfile is not an output from MuSICC grouping process. Please ensure you have the right file.\n";
    }
    $SVMclass{$gcells[0]}{'group'} = $gcells[4];
  }
  $line = $line + 1;
}

close GROUP;

#Write results to file
if(-e $output1)
{
  unlink $output1;
}
open WRITER, ">>$output1" or die "Cannot create $output1.\n";
my $writer = "SequenceId\tSVMClass\tProbability\tGroup\n";
print $writer;
print WRITER $writer;

foreach my $a (0 .. $#SequenceID)
{
  my $id = $SequenceID[$a];
  $writer = $id . "\t" . $SVMclass{$id}{'Class'} . "\t" . $SVMclass{$id}{'pval'} . "\t" . $SVMclass{$id}{'group'} . "\n";
  print $writer;
  print WRITER $writer;
}

close WRITER;

if(-e $output2)
{
  unlink $output2;
}
open WRITER, ">>$output2" or die "Cannot create $output2.\n";
$writer = "IonChan\tCount\n";
print "\n" . $writer;
print WRITER $writer;

foreach my $b (0 .. $#ionFam)
{
  my $ion = $ionFam[$b];
  $writer = $ion . "\t" . $ionCount{$ion} . "\n";
  print $writer;
  print WRITER $writer;
}

$writer = "\nTotal\t" . $totalChan . "\n";
print $writer;
print WRITER $writer;

close WRITER;

sub message {
  print "\nCompileResults.pl @ Bahiyah Nor, Sept 15th 2015\n\n";
  print " Usage: CompileResults.pl -s <SVMClassification.res> -g <sequencegrouping>\n";
  print " Options: -o <Output.tag> All the outputs from this program will be tagged as such\n";
  print "          -h For help\n";
  print " Compile the results of MuSICC run\n";
  print " Organise the classification made through the SVM classifier\n";
  print " with the groupings made\n";
  print " Present the results in a table with sequence id, classification, probability value\n";
  print " and the grouping (refer to manual regarding grouping)\n";
  die " Outputs: 1. MuSICC.[launchtime].Classification.txt 2. MuSICC.[launchtime].ChannelCount.txt\n\n";
  
}
