
MuSICC @ Multi-Screening Ion Channel Classifier
-- Bahiyah Nor, August 2015

** REQUIREMENTS
-------------------------------------
MuSICC is built using Perl scripts for Unix platform only.
The following softwares/packages/database are required for this program to run
  1. BLASTp
  2. InterProScan version 5.7.48
  3. TMHMM 2.0
  4. KEGG database
  5. R version 3.2.0 with 'e1071' package installed
We have only tested the program with the stated software versions.
We cannot guarantee that it will run perfectly with other versions.


** HOW TO RUN MuSICC
-------------------------------------
1.  BLASTp against KEGG database
    == $ blastp -query input_sequence -db KEGG_database -evalue 1e-15 -outfmt 6 -out blastp_output
    !! MuSICC requires tab-delimited output file, it will not process other output format
    
2.  Filter the BLASTp results
    == $ perl FirstFilter.pl -i blastp_output
    !! Two file outputs: 1. Headers of the filtered sequences 2. The report
    
3.  Retrieve the filtered sequences
    == $ perl RetrieveSeqs.pl -f original_fasta Output/FirstFilter.Headers.txt > FirstScreening.fa
    !! FirstFilter.Headers.txt is a file output from FirstFilter.pl. It is usually in Output/ folder.
    
4.  BLASTp against MuSICC.IonChannels.fa
    == $ blastp -query FirstScreening.fa -db Databases/MuSICC.IonChannels.fa -evale 1e-45 -outfmt 6 -out output_file
    !! MuSICC requires tab-delimited output file, it will not process other output format

5.  Filter the BLASTp results
    == $ perl SecondFilter.pl -i second_blastp_output
    !! A report of this process is produced

6.  Retrieve the filtered sequences
    == $ awk '{print $1}' second_blastp_output | uniq > Second_Headers.txt
       $ perl RetrieveSeqs.pl -f FirstScreening.fa Second_Headers.txt > FilteredSequences.fa

7.  Conserved Domain Predictions
    == $ interproscan.sh -i FilteredSequences.fa -f tsv --goterms -dp -o ipro.tsv
    !! interproscan.sh is a shell script that comes with InterProScan

8.  Transmembrane Domain Predictions
    == $ tmhmm -short FilteredSequences.fa > tmhmm.report.txt

9.  Sequence Grouping
    == $ perl SequenceGrouping.pl -i SecondFilter.report.txt -c ipro.tsv -i tmhmm.report.txt
    !! Output file: Output/SequenceGrouping.Report.txt

10. Classification via SVM classifier
    == $ perl PseAAC.DataPreProcess.pl -i FilteredSequences.fa -o FormattedFasta.fa
       $ perl BuildPseAAC.pl -i FormattedFasta.fa
       $ Rscript SVMClassifier.r

11. Produce the report %%Optional
    == $ perl CompileResults.pl -s SVMclassification_results -g SequenceGrouping.Report.txt

-- MuSICC allows BLASTp, InterProScan and TMHMM to be run independent of the pipeline
-- (i.e. it can all be run prior to the screening processes based on the initial fasta file)


** FAQ
-------------------------------------
1. Q: I have an initial BLASTp result, do I still need to run BLASTp?
   A: You don't have to run the first BLASTp (against the KEGG database) if you wish not to.
      But you should blast your sequences against MuSICC ion channel database.
2. Q: Is it necessary to blast my sequences against KEGG? Can I blast it against other database?
   A: Yes, you can. However, FirstFilter.pl cannot help you to filter your sequences as it is
      dependent on KEGG database sequence identifiers.
3. Q: I have a different version of InterProScan. Can I use that instead?
   A: As long as the output is in the same format, we don't think that this will cause any problem.
4. Q: Can I run MuSICC on different cores?
   A: No. But you can choose to run BLASTp on as many threads.
   

** ABOUT
-------------------------------------
MuSICC is a multi-screening ion channel classification pipeline.
It combines 3 publicly available bioinformatic tools:
  1. BLASTp from the BLAST+ suite
  2. InterProScan version 5.7.48
  3. TMHMM 2.0
The pipeline is divided into two parts, identification and classification.
Identification of ion channel sequences is done based on results obtained 
from the listed bioinformatic tools. Classification in ion channel families
is done by the built-in SVM classifier.

