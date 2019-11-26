# Base editor design tool
<p>This code designs every possible guide in the CDS region, 30 nucleotides into the intronic and UTR regions for 
user-defined transcripts. It then annotates the possible edits for each guide. Separate files annotating ClinVar SNPs are also generated.</p>
<b>Author</b>: Mudra Hegde, Ruth Hanna <br/>
<b>Email</b>: mhegde@broadinstitute.org, rhanna@broadinstitute.org <br/>
<b>Version: 2.0 </b> 

<b>Inputs</b>
1. <b>Input File</b>:.txt file with list of Ensembl transcript IDs in the first column and gene symbols in the second column OR 
FASTA file with nucleotide sequence
2. <b>Variant File</b>: Variant file: File with ClinVar SNPs (variant_summary.txt). This file can be downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz.
3. <b>Input type</b>: Indicate whether the file contains a list of transcripts or a nucleotide sequence.
4. <b>Base editor type</b>: Indicate the type of base editor (Rees et al.,2018) for which designs are required. This choice dictates the choice of PAM, sgRNA length, editing window and type of edit.
5. <b>PAM</b>: PAM preference if BE type has not been selected; Default: NGG.
6. <b>Edit window</b>: Editing window relative to nucleotide position in sgRNA, if BE type has not been selected; Default: 4-8.
7. <b>sgRNA length</b>: Length of sgRNA excluding PAM sequence, if BE type has not been selected; Default:20.
8. <b>Edit</b>: Type of edit made by base editor, if BE type has not been selected; Default: all, annotates for both C->T and A->G edit.
9. <b>Intron buffer</b>: Number of bp into the intron to consider for guide design.
10. <b>Filter GC</b>: Whether to filter out edits in a GC motif.
11. <b>Output name</b>: Name for output folder.

