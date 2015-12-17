#SequenceGrouping.pl @ Bahiyah Nor, Aug 17th 2015
#
# Usage: SequenceGrouping.pl 
#	-i <SecondFilter.Report> 	  Plain text file containing the sequence IDs and BLAST hit
#	-c <PredictedConservedDomain>	  Results obtained from InterProScan, with Pfam codes
#	-t <PredictedTransmembraneDomain> Table output from TMHMM, with transmembrane domain counts
# 	Optional: -h	Help information
# The program requires all three information to run. Without [-i][-c][-t] the program will terminate.
# Please ensure that the required information are provided.
use strict;
use warnings;

my $headerfile = "";
my $cdfile = "";
my $tmfile = "";
my $matchARG = 0;

foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i') #The sequence ids
  {
    if($i+1 < @ARGV)
    {
      $headerfile = $ARGV[$i+1];
      $matchARG = $matchARG + 1;
    }
  }
  elsif($ARGV[$i] eq '-c') #The Conserved domains file
  {
    if($i+1 < @ARGV)
    {
      $cdfile = $ARGV[$i+1];
      $matchARG = $matchARG + 1;
    }
  }
  elsif($ARGV[$i] eq '-t') #The transmembrane domain file
  {
    if($i+1 < @ARGV)
    {
      $tmfile = $ARGV[$i+1];
      $matchARG = $matchARG + 1;
    }
  }
  elsif($ARGV[$i] eq '-h')
  {
    print "#SequenceGrouping.pl @ Bahiyah Nor, Aug 17th 2015\n";
    print "#\n";
    print "# Usage: SequenceGrouping.pl\n";
    print "#\t -i <SecondFilter.Report>\t\t Plain text file containing the sequence IDs and BLAST hit\n";
    print "#\t -c <PredictedConservedDomain>\t\t Results obtained from InterProScan, with Pfam codes\n";
    print "#\t -t <PredictedTransmembraneDomain>\t Table output from TMHMM, with transmembrane domain counts\n";
    print "#\t Optional: -h    Help information\n";
    print "# The program requires all three information to run. Without [-i][-c][-t] the program will terminate.\n";
    print "# Please ensure that the required information are provided.\n";
    die "\n";
  }
}

if(@ARGV == 0 or $matchARG < 3)
{
  print "#SequenceGrouping.pl @ Bahiyah Nor, Aug 17th 2015\n";
  print "#\n";
  print "# Usage: SequenceGrouping.pl\n";
  print "#\t -i <SecondFilter.Report>\t\t Plain text file containing the sequence IDs and BLAST hit\n";
  print "#\t -c <PredictedConservedDomain>\t\t Results obtained from InterProScan, with Pfam codes\n";
  print "#\t -t <PredictedTransmembraneDomain>\t Table output from TMHMM, with transmembrane domain counts\n";
  print "#\t Optional: -h    Help information\n";
  print "# The program requires all three information to run. Without [-i][-c][-t] the program will terminate.\n";
  print "# Please ensure that the required information are provided.\n";
  die "\n";
}

#Load the Conserved domain profiles
my $cdprofile = "Database/MuSICC.ConservedDomProfiles.txt";
open CD, $cdprofile or die "Cannot load $cdprofile. Please ensure that the file exists.\n";
my %CDPRO;
my $cd;

print "Loading the conserved domain profiles...\n";
while($cd=<CD>)
{
  chomp($cd);
  my ($pcode, $ionFam, $prob) = split(/\s+/,$cd);
  $CDPRO{$ionFam}{$pcode} = $prob;
}

close CD;

#Load the transmembrane domain profiles
my $tmprofile = "Database/MuSICC.TransDomProfiles.txt";
open TM, $tmprofile or die "Cannot load $tmprofile. Plead ensure that the file exists.\n";
my %TMPRO;
my $tm;
my %CellTMCount;

print "loading the transmembrane domain profiles...\n";
while($tm=<TM>)
{
  chomp($tm);
  my @cells = split(/\s+/,$tm);
  if($cells[0] eq 'IonFam') #The header of the transmembrane domain profiles
  {
    foreach my $c (1 .. $#cells)
    {
      $CellTMCount{$c} = $cells[$c];  #Save the domain count respective to the column index
    }
  }
  else
  {
    foreach my $x (1 .. $#cells)
    {
      my $tmcount = $CellTMCount{$x};
      $TMPRO{$cells[0]}{$tmcount} = $cells[$x]; #The probability of the ion channel family to have this domain count
    }
  }
}

close TM;

my %Groupings;

#Load the sequence ids
open SID, $headerfile or die "Cannot load $headerfile. Please ensure that the file exists.\n";
my $sid;

while($sid=<SID>)
{
  chomp($sid);
  my @frags = split(/\s+/,$sid);
  if(!exists $Groupings{$frags[0]})  #Ensures that the sequence is only accounted for once
  {
    $Groupings{$frags[0]}{'BLAST'} = $frags[1];
  }
}

close SID;

#Load the conserved domain prediction result
open INTERPRO, $cdfile or die "Cannot load $cdfile. Please ensure that the file exists.\n";
my $interpro;
my %SeqConservedDomains;

while($interpro=<INTERPRO>)
{
  chomp($interpro);
  my @columns = split(/\s+/,$interpro);
  if(exists $Groupings{$columns[0]})  #if the sequence is included post BLAST
  {
    if($columns[3] eq 'Pfam')   #The line contains the pfam code
    {
      $SeqConservedDomains{$columns[0]}{$columns[4]} = 1; #$SeqConservedDomains{seqId}{pfamCode} = 1
    }
  }
}

close INTERPRO;

#Load the TMHMM results
open TMHMM, $tmfile or die "Cannot load $tmfile\n";
my $tmhmm;
my %SeqTransDomains;

while($tmhmm=<TMHMM>)
{
  chomp($tmhmm);
  my @splitted = split(/\s+/,$tmhmm);
  if(exists $Groupings{$splitted[0]})
  {
    my $tmco = substr($splitted[4],8);
    $SeqTransDomains{$splitted[0]} = $tmco; #$SeqTransDomains{seqId} = count
  }
}

close TMHMM;

#Match the profiles
#For each sequence in the set, find the best matching conserved domain profile
#Determine if the transmembrane domain is within the range of predicted class based on conserved domains and BLAST
#Write the results to output file: SequenceGrouping.Report.txt
my $output = "Output/SequenceGrouping.Report.txt";
if(-e $output)
{
  unlink $output;
}
open WRITER, ">>$output" or die "Cannot create $output.\n";
my $writer = "SequenceID\tBLAST_Hit\tCD_MatchedPro\tTM_MatchedPro\tGroup\n";
print WRITER $writer;

foreach my $seqId (keys %Groupings)
{
  print "Grouping $seqId...\n";

  #Find the Ion Channel subfamily that best matched the sequence  
  my $maxMatch = 0;
  my $CDclass = "NA";
  foreach my $family (keys %CDPRO)
  {
    my $totalMatch = 0;
    foreach my $pfam (keys %{ $SeqConservedDomains{$seqId} })
    {
      if(exists $CDPRO{$family}{$pfam})
      {
        if($CDPRO{$family}{$pfam} > 0)
        {
          $totalMatch = $totalMatch + 1;
        }
      }
    }
    #Determine if the sequence contains the common conserved domain in the class
    foreach my $dom (keys %{ $CDPRO{$family} })
    {
      if($CDPRO{$family}{$dom} > 0.95) #Common domain threshold = 0.95
      {
        $totalMatch = 0; #Reset the number of match to 0 if it is missing this domain
      }
    }
    
    if($totalMatch > $maxMatch)
    {
      $maxMatch = $totalMatch;
      $CDclass = $family;
    }
  }
  $Groupings{$seqId}{'ConDom'} = $CDclass;

  #Determine if the transmembrane domain count is within the range
  #Base on the subfamilies predicted by BLAST/Conserved Domains
  my $TMclass = "NA";
  my $blastPred = $Groupings{$seqId}{'BLAST'};
  my $seqTmCo = $SeqTransDomains{$seqId};
  if(exists $TMPRO{$CDclass}{$seqTmCo})
  {
    if($TMPRO{$CDclass}{$seqTmCo} > 0.05) #Common transmembrane domain count threshold = 0.05
    {
      $TMclass = "Matched";
    } 
  }
  if(exists $TMPRO{$blastPred}{$seqTmCo})
  {
    if($TMPRO{$blastPred}{$seqTmCo} > 0.05)
    {
      $TMclass = "Matched";
    }
  }
  $Groupings{$seqId}{'TransDom'} = $TMclass;

  #Group the sequence
  if($CDclass ne 'NA')
  {
    if($TMclass ne 'NA')
    {
      $Groupings{$seqId}{'Group'} = 1;
    }
    else
    {
      $Groupings{$seqId}{'Group'} = 2;
    }
  }
  else
  {
    if($TMclass ne 'NA')
    {
      $Groupings{$seqId}{'Group'} = 3;
    }
    else
    {
      $Groupings{$seqId}{'Group'} = 4;
    }
  }

  #Write to file
  $writer = $seqId . "\t" . $Groupings{$seqId}{'BLAST'} . "\t" . $CDclass . "\t" . $TMclass . "\t" . $Groupings{$seqId}{'Group'} . "\n";
  print WRITER $writer;

}

close WRITER;
