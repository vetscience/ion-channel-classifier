#BuildPseAAC.pl @ Bahiyah Nor, May 22 2015
#
# Usage: BuildPseAAC.pl -i <fastafile> 
# Other option: -l <lambda> Lambda value to be used according to Chou's PseAAC, [l > 0]
#               -w <weight> The weight to be used according to Chou's PseAAC,[0.05 < w < 0.7]
#               -h          For help information
# 
# Build the PseAAC of the data provided based on Chou's PseAAC (REF) and the dipeptide composition
# using the hydropohibicity, hydrophilicity and mass of the side chain
# assuming the the file is in a fasta format
# **Lambda value: has to be a non-negative integer and lesser than the length of the shortest sequence
# **Weight: has to be between 0.05 and 0.7
use strict;
use warnings;

my $input = "";
my $output = "Output/PseAAC.csv";
if(-e $output)
{
  unlink $output;
}
my $lambda = 55; #default lambda, the program may change the value automatically depending on the length of sequences
my $weight = 0.7; #default weight
my $matchArg = 0;
my $inputFlag = 0;
my $shortest = 100000000;

foreach my $i (0 .. $#ARGV)
{
  if($ARGV[$i] eq '-i')
  {
    if($i+1 < @ARGV)
    {
      $input = $ARGV[$i+1];
      open TMP, $input or die "Cannot load $input.\n";
      #Adjust the lambda value depending on the length of the shortest sequence in the set
      my $line;
      while($line=<TMP>)
      {
        chomp($line);
        if(substr($line,0,1) ne '>')
        {
          if(length($line) < $shortest)
          {
            $shortest = length($line);
          }
        }
      }
      if($shortest < $lambda)
      {
        $lambda = $shortest - 2;
        print "Lambda reset to $lambda\n";
      }
      close TMP;
      $matchArg = 1;
      $inputFlag = 1;
    }
    else
    {
      print "\n**!! Error encountered. The input file is not specified.\n\n";
    }
  }
  elsif($ARGV[$i] eq '-l')
  {
    if($i+1 < @ARGV)
    {
      if($ARGV[$i+1] < $shortest)
      {
        $lambda = $ARGV[$i+1];
        print "Lambda reset to $lambda\n";
      }
      else
      {
        print "The lambda value defined is larger than the shortest sequence in the set.\nThe default value will be used.\n";
      }
      $matchArg = 1;
    }
    else
    {
      print "Invalid lambda value. The default value will be used.\n";
    }
  }
  elsif($ARGV[$i] eq '-w')
  {
    if($i+1 < @ARGV)
    {
      if($ARGV[$i+1] < 0.7 and $ARGV[$i+1] > 0.05)
      {
        $weight = $ARGV[$i+1];
      }
      else
      {
        print "Invalid weight value. The default value will be used.\n";
      }
      $matchArg = 1;
    }
    else
    {
      print "Invalid weight value. The default value will be used.\n";
    }
  }
}

if(@ARGV == 0 or $inputFlag == 0 or $matchArg == 0)
{
  print "#BuildPseAAC.pl @ Bahiyah Nor, May 22 2015\n#\n";
  print "# Usage: BuildPseAAC.pl -i <fastafile>\n"; 
  print "# Other option: -l <lambda> Lambda value to be used according to Chou's PseAAC, [l > 0]\n";
  print "#               -w <weight> The weight to be used according to Chou's PseAAC,[0.05 < w < 0.7]\n";
  print "#               -h          For help information\n#\n"; 
  print "# Build the PseAAC of the data provided based on Chou's PseAAC and the dipeptide composition\n";
  print "# using the hydropohibicity, hydrophilicity and mass of the side chain\n";
  print "# assuming the the file is in a fasta format\n";
  print "# **Lambda value: has to be a non-negative integer and lesser than the length of the shortest sequence\n";
  die "# **Weight: has to be between 0.05 and 0.7\n\n";
}

#Declare the hydrophobicity, hydrophilicity and mass
my %composition;
my $psefile = "Database/PseAAC.table.txt";
open PSE, $psefile or die "Cannot load $psefile\n";
my $ps;
while($ps = <PSE>)
{
  chomp($ps);
  my ($aa, $compose) = split(/\t/, $ps, 2);
  $composition{$aa} = $compose;
}
#End of loading composition
close PSE;

#Convert the values into the standard version
my $HB_totdiv20 = 0;
my $HL_totdiv20 = 0;
my $M_totdiv20 = 0;
for(my $b=0; $b<=2; $b++)
{
  my $tot=0;
  foreach (keys %composition)
  {
    my $tmp = $composition{$_};
    my @all = split(/\t/, $tmp);
    $tot = $tot + $all[$b];
  }
  $tot = $tot / 20;
  
  my $denom = 0;
  foreach (keys %composition)
  {
    my $tmp2 = $composition{$_};
    my @all2 = split(/\t/, $tmp2);
    $denom = ($all2[$b] - $tot)^2;
  }
  $denom = sqrt($denom / 20);
  
  foreach(keys %composition)
  {
    my $tmp3 = $composition{$_};
    my @values = split(/\t/,$tmp3);
    my $val = $values[$b];
    $values[$b] = ($val - $tot) / $denom;
    my $back = $values[0] . "\t" . $values[1] . "\t" . $values[2];
    $composition{$_} = $back;
  }
}

#----------------------------------------------------------------------------------------------------------
#Chou's Pseudo Amino Acid Composition: 
#Compute the amino acid composition for the sequence
#using the formula:
#L = length(seq)
#P = [t1,t2,t2,...,t(lambda)]  where lambda < L;
#t1 = (1/(L-1)) * [sum of correlation function between r(i) and r(i+1)]
#t2 = (1/(L-2)) * [sum of correlation function between r(i) and r(i+2)]
#...
#t(lambda) = (1/(L-lambda))) * [sum of correlation function between r(i) and r(i+lambda)];
#
#Correlation function:
#CR(r[i],r[j]) = (1/3)*{(HB(r[j])-HB(r[i]))^2 + (HL(r[j])-HL(R[i]))^2 + (M(r[j]-M(r[i]))^2}
#   where HB = hydrophobicity, HL = hydropholicity, M = side-chain mass
#
#All the values were converted using:
# C1 = (C0 - {sum of (Ci/20)}) / sqrt({sum of (C0 - sum of[Ci/20])}/20)
#----------------------------------------------------------------------------------------------------------
open SEQ, $input or die "Cannot load $input\n";
my $read;
my %PseAAC;
for(my $k=1; $k<26; $k++)
{
  $PseAAC{$k} = -1;
}
my $prev_id = "";
my $firstline = 0;
my $header = "Gene	S	F	T	N	K	Y	E	V	Q	M	C	L	A	W	P	H	D	I	R	G	SS	SF	ST	SN	SK	SY	SE	SV	SQ	SM	SC	SL	SA	SW	SP	SH	SD	SI	SR	SG	FS	FF	FT	FN	FK	FY	FE	FV	FQ	FM	FC	FL	FA	FW	FP	FH	FD	FI	FR	FG	TS	TF	TT	TN	TK	TY	TE	TV	TQ	TM	TC	TL	TA	TW	TP	TH	TD	TI	TR	TG	NS	NF	NT	NN	NK	NY	NE	NV	NQ	NM	NC	NL	NA	NW	NP	NH	ND	NI	NR	NG	KS	KF	KT	KN	KK	KY	KE	KV	KQ	KM	KC	KL	KA	KW	KP	KH	KD	KI	KR	KG	YS	YF	YT	YN	YK	YY	YE	YV	YQ	YM	YC	YL	YA	YW	YP	YH	YD	YI	YR	YG	ES	EF	ET	EN	EK	EY	EE	EV	EQ	EM	EC	EL	EA	EW	EP	EH	ED	EI	ER	EG	VS	VF	VT	VN	VK	VY	VE	VV	VQ	VM	VC	VL	VA	VW	VP	VH	VD	VI	VR	VG	QS	QF	QT	QN	QK	QY	QE	QV	QQ	QM	QC	QL	QA	QW	QP	QH	QD	QI	QR	QG	MS	MF	MT	MN	MK	MY	ME	MV	MQ	MM	MC	ML	MA	MW	MP	MH	MD	MI	MR	MG	CS	CF	CT	CN	CK	CY	CE	CV	CQ	CM	CC	CL	CA	CW	CP	CH	CD	CI	CR	CG	LS	LF	LT	LN	LK	LY	LE	LV	LQ	LM	LC	LL	LA	LW	LP	LH	LD	LI	LR	LG	AS	AF	AT	AN	AK	AY	AE	AV	AQ	AM	AC	AL	AA	AW	AP	AH	AD	AI	AR	AG	WS	WF	WT	WN	WK	WY	WE	WV	WQ	WM	WC	WL	WA	WW	WP	WH	WD	WI	WR	WG	PS	PF	PT	PN	PK	PY	PE	PV	PQ	PM	PC	PL	PA	PW	PP	PH	PD	PI	PR	PG	HS	HF	HT	HN	HK	HY	HE	HV	HQ	HM	HC	HL	HA	HW	HP	HH	HD	HI	HR	HG	DS	DF	DT	DN	DK	DY	DE	DV	DQ	DM	DC	DL	DA	DW	DP	DH	DD	DI	DR	DG	IS	IF	IT	IN	IK	IY	IE	IV	IQ	IM	IC	IL	IA	IW	IP	IH	ID	II	IR	IG	RS	RF	RT	RN	RK	RY	RE	RV	RQ	RM	RC	RL	RA	RW	RP	RH	RD	RI	RR	RG	GS	GF	GT	GN	GK	GY	GE	GV	GQ	GM	GC	GL	GA	GW	GP	GH	GD	GI	GR	GG	L.1	L.2	L.3	L.4	L.5	L.6	L.7	L.8	L.9	L.10	L.11	L.12	L.13	L.14	L.15	L.16	L.17	L.18	L.19	L.20	L.21	L.22	L.23	L.24	L.25	L.26	L.27	L.28	L.29	L.30	L.31	L.32	L.33	L.34	L.35	L.36	L.37	L.38	L.39	L.40	L.41	L.42	L.43	L.44	L.45	L.46	L.47	L.48	L.49	L.50	L.51	L.52	L.53	L.54	L.55";

#my @splitHead = split(/ /,$header);
#$header = join("\t",@splitHead);

open WRITER, ">>$output" or die "Cannot create $output\n";
print WRITER "$header\n";
my $unknownflag = 0;

print "Computing lambda = $lambda\n";

while($read = <SEQ>)
{
  chomp($read);
  if(substr($read,0,1) eq '>')  #The header
  {
    my @words = split(' ', $read);
    if(($words[0] ne $prev_id) and ($firstline > 0) and ($unknownflag ==0))
    {
      my $write = $prev_id;
      

      foreach (keys %PseAAC)
      {
        $write = $write . "\t" . $PseAAC{$_};
      }
      $write = $write . "\n";
      print WRITER $write;
    }
    $firstline = $firstline + 1;
    $prev_id = $words[0];
    $unknownflag = 0;
  }
  elsif (length($read) > 1)
  {    
    my $seq_length = length($read);
    
    #Compute the frequency of the amino acids
    my %frequency;
    my $total_frequency = 0;
    my $aac = 1;    
    foreach (keys %composition)
    {
      my $cur_aa = $_;
      my $aa_count = 0;
      for(my $a=0; $a<$seq_length; $a++)
      {
        if(substr($read,$a,1) eq $cur_aa)
        {
          $aa_count = $aa_count + 1;
        }
      }
      $frequency{$aac} = $aa_count / 20;
      $total_frequency = $total_frequency + $frequency{$aac};
      $aac = $aac+1;
    }
    
    #Compute the dipeptide composition !!!There are 400 of these
    #AB =/= BA
    my %dipeptide;
    my $totDiPeptideFreq = 0;
    my $totNumDP = $seq_length - 1;
    my $dp = 1;
    foreach my $first (keys %composition)
    {
      foreach my $second (keys %composition)
      {
        my $di = $first . $second;
        my $di_count = 0;
      
        for(my $d=0; $d<$seq_length-1; $d++)
        {
          if(substr($read,$d,2) eq $di)
          {
            $di_count = $di_count + 1;
          }
        }
        $totDiPeptideFreq = $totDiPeptideFreq + $di_count;
        $dipeptide{$dp} = $di_count / $totNumDP;
        $dp = $dp + 1;
      }
    }
    
    #compute each theta value
    my %theta;
    my $total_theta = 0;
    for(my $t=1; $t<=$lambda; $t++)
    {
      my $total_sum = 0;
            
      for(my $i=0; $i<($seq_length-$t-1); $i++)
      {        
        my $Ri = substr($read,$i,1);
        my $Rj = substr($read,$i+$t,1);
        
        my @ri_compose;
        my @rj_compose;
        my $tmpri;
        if(exists $composition{$Ri})
        {
          $tmpri = $composition{$Ri};
          @ri_compose = split(/\t/, $tmpri);
        }
        else
        {
          #If the alphabet is X, we will compute using A composition (arbitarily chosen)
          if($Ri eq 'X')
    	    {
	          $tmpri = $composition{'A'};
	        }
          elsif($Ri eq 'B')
	        {
	          $tmpri = $composition{'D'};
          }
          elsif($Ri eq 'Z')
          {
	          $tmpri = $composition{'Q'};
          }
	        else
          {
            $unknownflag = 1;
          }
          if($unknownflag == 0)
 	        {
            @ri_compose = split(/\t/,$tmpri); 
	        }
        }
        
        my $tmprj;
        if(exists $composition{$Rj})
        {
          $tmprj = $composition{$Rj};
          @rj_compose = split(/\t/, $tmprj);
        }
        else
        {
	        if($Rj eq 'X')
	        {
	          $tmprj = $composition{'A'};
	        }
	        elsif($Rj eq 'B')
	        {
	          $tmprj = $composition{'D'};
	        }
	        elsif($Rj eq 'Z')
	        {
	          $tmprj = $composition{'Q'};
	        }
	        else
	        {
            $unknownflag = 1;
      	  }
	        if($unknownflag == 0)
	        {
	          @rj_compose = split(/\t/, $tmprj);
	        }
        }
        
	      if($unknownflag == 0)
        {
          my $correlation = 0;
          for(my $x=0; $x < 3; $x++)
          {
            $correlation = $correlation + (($rj_compose[$x] - $ri_compose[$x])^2);
          }
        
          my $tt = (1/3) * $correlation;
          $total_sum = $total_sum + $tt;
 	      }
      }
      
      if($unknownflag == 0)
      {
 	      $theta{$t} = (1/($seq_length-$t)) * $total_sum;
        $total_theta = $total_theta + $theta{$t};
      }
    }# End of computing theta
    
    if($unknownflag == 0)
    {
      #Compute the PseAAC
      my $bottom = $total_frequency + ($weight*$total_theta);
      for(my $z=1; $z<=(420+$lambda); $z++)
      {
        my $tmps;
        if($z <= 20)
        {
          $tmps = $frequency{$z} / $bottom;
        }
        elsif($z > 20 and $z <= 420)
        {
          $tmps = $dipeptide{$z - 20};
        }
        else
        {
          $tmps = $theta{$z-420} / $bottom;
        }
        $PseAAC{$z} = $tmps;
      }
    } 
  }
}
#end of reading SEQ

if($unknownflag == 0)
{
  my $like = $prev_id;
  foreach (keys %PseAAC)
  {
    $like = $like . "\t" . $PseAAC{$_};
  }
  $like = $like . "\n";
  print WRITER $like;
}

close SEQ;
