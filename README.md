ExRec:  A PYTHON PIPELINE FOR GENERATING RECOMBINATION-FILTERED MULTI-LOCUS DATASETS
Sam McCarthy Potter1 and W. Bryan Jennings2
1Department of Biology, Carleton College, Northfield, Minnesota 55057, USA.                                  Email: smccarthypotter@gmail.com
2 Department of Evolution, Ecology, & Organismic Biology, University of California, Riverside, California 92521, USA. Email: bryan.jennings@ucr.edu

Table of Contents
S1. About this software package									1
S2. Conventions followed in this user manual							2
S3. Installation procedures									2
S4. Overview of running the program scripts from the command line				2 
S5. Detailed instructions for running each program script from the command line			3
S6. Testing the program scripts using the example data						12
S7. Supplementary references									13

S1. About this software package
The ExRec package allows users to easily obtain recombination-filtered datasets that are often required in coalescent-based phylogeography, population genomics, and shallow-scale phylogenomics studies. The package consists of five stand-alone python scripts, which form a three-step pipeline (step 3 is optional). In step 1, the user executes a program script on the command line to batch-convert single-locus NEXUS or PHYLIP files into a concatenated partitioned interleaved NEXUS file, which is the standard input file for the package’s main program script. In step 2, the user inputs the concatenated NEXUS file into the main program, which performs the following steps in automatic fashion: 1) conducts four-gamete tests (Hudson and Kaplan, 1985) on all loci; 2) truncates apparently recombined loci down to non-recombined sequence blocks (following Hey and Nielsen, 2004); and 3) outputs recombination-filtered NEXUS and PHYLIP files and a tab-delimited summary table that shows descriptive statistics for each locus and the results of the recombination-filtering analyses. The user can select whether to output the longest (following Hey and Nielsen, 2004) or randomly selected (following Hey and Wang, 2019) non-recombined block for each locus that showed evidence of one or more historical recombination events. The concatenated-loci output files, which are in NEXUS and PHYLIP formats, can then be used in popular software programs like BPP (Yang, 2015; Flouri et al., 2018). In the optional step 3, two program scripts in the package allow the user to split the recombination-filtered concatenated data files into single-locus files in NEXUS and PHYLIP formats, which can be batch-input into phylogenetics programs thus facilitating summary methods species tree analyses (e.g., ASTRAL Mirarab et al., 2014, 2016).
S2. Conventions followed in this user manual
The term “program script” refers to each of the five stand-alone Python 3 programs included in the ExRec package: Nexcombine.py, Phycombine.py, FGT.py, Nexsplit.py, and Physplit.py. The names of these scripts are italicized throughout the document. The user runs each of the five scripts from the command line. All command line syntax instructions below are shaded in gray. Datafile names are underlined, and folder names are in bold letters.

S3. Installation procedures
Python 3 must be installed on your computer. If your computer does not have Python 3, then you can install it from https://www.python.org/downloads/
Download the ExRec package from https://github.com/Sammcarthypotter/ExRec. The package includes: five program scripts (i.e., Nexcombine.py, Phycombine.py, FGT.py, Nexsplit.py, and Physplit.py) and two test data sets (i.e., “finch loci” and “hominoid loci”). These programs can run on Macintosh, PC, or UNIX/LINUX computers.
Place the Nexcombine.py, Phycombine.py, FGT.py, Nexsplit.py, and Physplit.py program scripts into separate folders on your desktop or whichever location on your computer that you choose.

S4. Overview of running the program scripts from the command line
Place your single-locus NEXUS files into a folder named nexus_files or your single-locus PHYLIP files into a folder named phylip_files. If you are working with NEXUS files, then place the folder nexus_files into the folder that contains Nexcombine.py. If instead you are working with PHYLIP files, then place phylip_files into the folder that contains Phycombine.py. Execute Nexcombine.py or Phycombine.py by entering the appropriate commands on the command line prompt (explained in section S5). After Nexcombine.py or Phycombine.py finishes, an output file named combineloci.nex will appear in the folder containing the program script and data folder. This file contains all loci together in a concatenated partitioned interleaved NEXUS file, which is the required input file for FGT.py.
Next, move the combineloci.nex file to the folder containing the FGT.py script. Execute FGT.py by entering the appropriate commands on the command line prompt (explained in section S5). When FGT.py finishes it run, it will output three files: a recombination-filtered concatenated partitioned interleaved NEXUS file named trunc_combineloci.nex, a recombination-filtered concatenated PHYLIP file named trunc_combineloci.phy, and a tab-delimited text file named Results Summary Table.txt.
If you would like to split the trunc_combineloci.nex file into single-locus NEXUS files, then place the trunc_combineloci.nex file inside the same folder as Nexsplit.py. If instead you would like to split the data in trunc_combineloci.nex into single-locus PHYLIP files, then place the trunc_combineloci.nex inside the same folder as Physplit.py. Execute Nexsplit.py or Physplit.py by entering the appropriate commands on the command line prompt (explained in section S5). When Nexsplit.py finishes, it should output a folder named nexus_split_files that contains all single-locus NEXUS files and when Physplit.py finishes it should output a folder named phylip_split_files that contains all single-locus PHYLIP files.



S5. Detailed instructions for running each program script from the command line
1. The Nexcombine.py script
The Nexcombine.py script converts a collection of single-locus NEXUS files into one concatenated NEXUS file—the standard input file for FGT.py. To use Nexcombine.py, you must input your data as single-locus NEXUS files in either “sequential” or “interleaved” formats. Below is an example of a single-locus NEXUS sequential file (note: most of the sequences are not shown):


Example of a single-locus NEXUS interleaved file (bottom of file not shown):


Running Nexcombine.py on the command line
To run Nexcombine.py, you must first place the Nexcombine.py script and a folder containing the single-locus NEXUS files into the same folder. The folder containing the single-locus NEXUS files must be named nexus_files and there cannot be any other files present in the folder. Next, cd to the folder containing Nexcombine.py and nexus_files. On the command line, enter the following commands to run the program:
>python3 Nexcombine.py
After the program starts running, it will print the names of all the single-locus NEXUS files on the screen and then list the first and last sites of each locus in brackets, thus giving the length (in bp) of each locus. The output file will be named combineloci.nex, which contains the concatenated partitioned loci in NEXUS interleaved format. An example combineloci.nex file is shown below (note: only the top part of file is showing):


You will see that the bottom of the combineloci.nex file contains a “charset block,” which lists all loci by their names (within the single quotation marks) and their coordinates within the data matrix:


Before using your combineloci.nex file in FGT.py, it is very important to open this file in a text editor and check that it is free of any obvious errors. One simple check that can be done is to open the first two single-locus NEXUS files and compare them to the first two loci in the combineloci.nex file, ensuring that the new file contains the correct “ntaxa,” “nchar,” and “charset” values and shows full length sequences.


2. The Phycombine.py script
The Phycombine.py script converts a collection of single-locus PHYLIP files into a concatenated NEXUS file named combineloci.nex for input into FGT.py. Although all PHYLIP files show the number of sequences and the number of sites at the top of the files followed by the aligned sequences, it is important to understand that many variants of this format exist. PHYLIP files either have sequence names in the “strict” format (i.e., maximum of 10 characters long) or “relaxed” format (i.e., names can be up to 256 characters long) and they can be in either the “sequential” or “interleaved” formats. Accordingly, to ensure that Phycombine.py produces an error-free combineloci.nex file, you must use the correct command line syntax depending on which PHYLIP variant files are being input into Phycombine.py. Table 1 lists various PHYLIP file types that are compatible with Phycombine.py as well as the command line syntax required to process each one.

Table 1. Summary of PHYLIP format variants that Phycombine.py can process. Flags can be lower- or upper-case and their order of input on the command line does not affect the analysis. Note, other PHYLIP format variants may work in Phycombine.py but you should check output files for errors.
PHYLIP format variant
required flags
command line syntax
STRICT SEQUENTIAL
none (defaults)
>python3 Phycombine.py
STRICT SEQUENTIAL UPPER
u
>python3 Phycombine.py u
STRICT SEQUENTIAL WRAPPED
w
>python3 Phycombine.py w
STRICT SEQUENTIAL UPPER WRAPPED
u w
>python3 Phycombine.py u w
STRICT INTERLEAVED
none (defaults)
>python3 Phycombine.py
STRICT INTERLEAVED GAPPED
none (defaults)
>python3 Phycombine.py
RELAXED SEQUENTIAL
r
>python3 Phycombine.py r
RELAXED SEQUENTIAL WRAPPED
r w
>python3 Phycombine.py r w
RELAXED INTERLEAVED
r
>python3 Phycombine.py r
RELAXED INTERLEAVED GAPPED
r
>python3 Phycombine.py r



Descriptions of each PHYLIP format variant are provided below. Only the beginning part of each file type is shown.

STRICT SEQUENTIAL: sequence name is in PHYLIP strict name format followed by the sequence in sequential format. The sequence must start at the 11th character space. Required flags: none (defaults)



STRICT SEQUENTIAL UPPER: sequence name is in PHYLIP strict name format on one line and then the sequence starts at the beginning of the next line in sequential format. Required flags: u


STRICT SEQUENTIAL WRAPPED: sequence name is in PHYLIP strict name format followed by the sequence in sequential format. The sequence starts at the 11th character space and then wraps around to the lines below. Required flags: w







STRICT SEQUENTIAL UPPER WRAPPED: sequence name is in PHYLIP strict name format on one line and then the sequence starts at the beginning of the next line in sequential format. Sequence wraps around to the lines below. Required flags: u w


STRICT INTERLEAVED: sequence name is in PHYLIP strict name format followed by the sequence in interleaved format. The sequence must start at the 11th character space. Required flags: none (defaults)

STRICT INTERLEAVED GAPPED: sequence name is in PHYLIP strict name format followed by the sequence in interleaved format. The sequence must start at the 11th character space. The sequences contain a single-character gap after every ten nucleotides. The sequence on each line can be any length. Required flags: none (defaults)


RELAXED SEQUENTIAL: sequence name is in PHYLIP relaxed name format followed by at least one space followed by the sequence in sequential format. Required flags: r


RELAXED SEQUENTIAL WRAPPED: sequence name is in PHYLIP relaxed name format followed by at least one space followed by the sequence in sequential format. Sequence wraps around to the lines below. Required flags: r w

RELAXED INTERLEAVED: sequence name is in PHYLIP relaxed name format followed by at least one space followed by the sequence in interleaved format. Sequence on each line can be any length (e.g., 50 nucleotides long). Required flags: r


RELAXED INTERLEAVED GAPPED: sequence name is in PHYLIP relaxed name format followed by at least one space followed by the sequence in interleaved format. The sequences contain a single-character gap after every ten nucleotides. The sequence on each line can be any length. Required flags: r


Running Phycombine.py
To run Phycombine.py, you must first place the Phycombine.py script and a folder containing the single-locus PHYLIP files into the same folder. The folder containing the single-locus files must be named phylip_files and there cannot be any other files present in this folder. Next, cd to the folder containing Phycombine.py and phylip_files. On the command line, enter the commands to run the program using the syntax provided in Table 1 as a guide. After starting Phycombine.py, the program will first print the names of the single-locus PHYLIP files to the screen and then it will print the name of the program followed by “STRICT,” “STRICT UPPER,” “STRICT WRAPPED,” “STRICT UPPER WRAPPED,” RELAXED, or RELAXED WRAPPED and the flags used to indicate the program settings for the input files. When finished, Phycombine.py will output the data in the combineloci.nex file. Before using this newly generated combineloci.nex file in FGT.py, it is very important to first open the file in a text editor and check that it is free of any obvious errors. Incorrectly generated combineloci.nex files may lack most or all the sequences, or the sequences may be oddly arranged. You can also easily spot flawed combineloci.nex files if the numbers of sequences and/or sites for each locus are wrong. These flawed files, which in our experience are easy to notice, can be produced when we input the wrong command line syntax for the input files. Thus, it is essential that you use the correct command line syntax when running Phycombine.py.
3. The FGT.py script
The FGT.py script automatically carries out the four-gamete filtering procedures on all input loci in the combineloci.nex file and then outputs three files: a recombination-filtered concatenated partitioned interleaved NEXUS file named trunc_combineloci.nex, a recombination-filtered concatenated PHYLIP file named trunc_combineloci.phy, and a tab-delimited text file named Results Summary Table.txt. FGT.py can perform two different four-gamete filtering analyses: 1) it can output the longest presumably non-recombined sequence block for each apparently recombined locus (current best practice, first suggested by Hey and Nielsen, 2004); or 2) it can output a randomly selected non-recombined sequence block when there are at least two presumably non-recombined blocks for an apparently recombined locus (an approach first suggested by Hey and Wang, 2019).

Option 1: output the longest non-recombined sequence blocks (“longest block mode”)
To run FGT.py in the mode that outputs the longest non-recombined sequence blocks, you first need to place the combineloci.nex file in the same folder as the FGT.py script. Then cd to the folder containing FGT.py and combineloci.nex before entering the following commands on the command line:
>python3 FGT.py combineloci.nex ? – ms
After entering the program name “FGT.py” on the command line, you must then enter a single space followed by the name of the input file, which will usually be combineloci.nex. You can input other file names provided that the file name exactly matches the name of the file. After entering the name of the input file, you must then enter a single space, a question mark, another single space, a hyphen, another single space, and the letters “ms.” The “?”, “–”, and “ms” flags instruct FGT.py to carry out the four-gamete tests under the most stringent conditions because most datasets will have missing data at some sites (indicated by “?” in the matrices), alignment gaps (indicated by “–” in the matrices) and will be missing entire sequences for some loci (indicated by sequences of question marks). Thus, the “?” and “–” flags instruct FGT.py to ignore sites with missing data and alignment gaps (following Rozas et al., 2017), respectively, while the “ms” flag instructs the program to ignore missing sequences while conducting four-gamete tests. These settings should yield identical four-gamete test results to those obtained using the program DNAsp (Rozas et al., 2017).
While FGT.py is running, each locus name will appear on the computer screen in the order that it is processed. When FGT.py finishes, a timer on the screen will indicate how long the program ran and three files will be output into the FGT.py folder: trunc_combineloci.nex, trunc_combineloci.phy, and Results Summary Table.txt. The recombination-filtered data are thus contained in the concatenated NEXUS and PHYLIP files trunc_combineloci.nex and trunc_combineloci.phy, respectively. You can utilize Nexsplit.py or Physplit.py in the ExRec package to split the trunc_combineloci.nex file into its component single-locus NEXUS or PHYLIP files, respectively, for use in batch-mode phylogenetic analyses thereby facilitating summary methods species tree analyses (see below). You can also immediately input the trunc_combineloci.phy file into species delimitation/historical demography analyses (e.g., the software BPP; Yang, 2015; Flouri et al., 2018).
Results Summary Table.txt is a text file that contains descriptive statistics about each locus and the results of the four-gamete tests/recombination filtering procedures. The table shows for each locus: locus name, starting length (bp), length (bp) excluding gaps/missing data at sites, S (number of segregating sites), list of sites that violate the infinite sites model, RM (minimum number of recombination events; Hudson and Kaplan, 1985), pairs of sites that had recombination event(s) within them, sites that define the longest non-recombined block, and the length of the retained longest non-recombined block (bp). You can copy and paste this file into a spreadsheet for further metanalyses of the data (Fig. 1).

Fig. 1. Summary table of descriptive statistics for all loci and results when FGT.py is run in longest block mode.

Option 2: output randomly chosen non-recombined sequence blocks (“random block mode”)
To run FGT.py in the mode that outputs randomly chosen non-recombined sequence blocks, you first need to place the combineloci.nex file in the same folder as the FGT.py script. Then cd to the folder containing FGT.py and combineloci.nex before entering the following commands on the command line:
>python3 FGT.py combineloci.nex ? – ms r
Note that you must enter an “r” flag on the command line after you input the file name and any flags that are used. As before, single spaces must separate each flag. When you run FGT.py in the random block mode, the summary table will be identical to the table generated under the longest block mode except for the following differences: columns showing the pairs of sites that define presumably non-recombined sequence blocks, the pair of sites that define the randomly selected non-recombined block, and the length (bp) of the randomly selected block (Fig. 2).

Fig. 2. Summary table of descriptive statistics for all loci and results when FGT.py is run in random block mode.

4. The Nexsplit.py script
The Nexsplit.py script converts the trunc_combineloci.nex file into single-locus NEXUS files in interleaved format. To use Nexsplit.py, you must place the trunc_combineloci.nex file into the same folder as the Nexsplit.py script. After you cd to the folder containing Nexsplit.py and trunc_combineloci.nex, enter the following commands on the command line:
>python3 Nexsplit.py
While the program is running, it will print the name of the input file in brackets on the screen. When finished, Nexsplit.py will output single-locus NEXUS files into a new folder named nexus_split_files.

5. The Physplit.py script
The physplit.py script converts the trunc_combineloci.nex file into single-locus PHYLIP files in relaxed sequential format. To use Physplit.py, you must place the trunc_combineloci.nex file into the same folder as the Physplit.py script. After you cd to the folder containing Physplit.py and trunc_combineloci.nex, enter the following commands on the command line:
>python3 Physplit.py
While the program is running, it will print the name of the input file in brackets on the screen. When finished Physplit.py will output single-locus PHYLIP files into a new folder named phylip_split_files.

S6. Testing the program scripts using the example data
You can test the five program scripts using two example datasets that are also provided with the ExRec package.
The finch_nexus folder contain 27 anonymous DNA sequence loci from the study of Jennings and Edwards (2005). The single-locus files are in NEXUS sequential format and thus should be input into Nexcombine.py using the above instructions. You will need to rename this folder nexus_files so it can be used in Nexcombine.py.
The hominoids_phylip folder contains 292 anonymous DNA sequence loci from the study of Costa et al. (2016). The single-locus files are in PHYLIP strict sequential upper format and thus should be input into Phycombine.py. Remember to rename the input data folder as phylip_files before running the program.
When finished running, Nexcombine.py and Phycombine.py should each output a single file that contains all loci in a concatenated partitioned interleaved NEXUS format named combineloci.nex. You can then immediately use this file as the input data file for FGT.py using the above instructions.
When FGT.py finishes its run, it should output three files: a recombination-filtered concatenated partitioned interleaved NEXUS file named trunc_combineloci.nex, a recombination-filtered concatenated PHYLIP file named trun_combineloci.phy, and a tab-delimited text file named Results Summary Table.txt that summarizes the loci characteristics and analysis results. You can copy and paste the text file into a spreadsheet for easier visualization of the metadata.
You can then test Nexsplit.py by inputting the trunc_combineloci.nex file into Nexsplit.py. When finished the program should output the data into a folder named nexus_split_files, which contains single-locus sequential NEXUS files.
Finally, you can test your copy of Physplit.py by inputting the trunc_combineloci.nex file into Physplit.py. When finished running, this program should output a folder named phylip_split_loci, which contains the single-locus PHYLIP files.

S7. Supplementary references
Costa,I.R. et al. (2016) In silico phylogenomics using complete genomes: a case study on the evolution of hominoids. Genome Res., 26, 1257–1267.
Flouri,T. et al. (2018) Species tree inference with BPP using genomic sequences and the multispecies coalescent. Molecular Biology and Evolution, 35, 2585–2593.
Hey,J. and Nielsen,R. (2004) Multilocus methods for estimating population sizes, migration rates and divergence time, with applications to the divergence of Drosophila pseudoobscura and D. persimilis. Genetics, 167, 747–760.
Hey,J. and Wang,K. (2019) The effect of undetected recombination on genealogy sampling and inference under an isolation‐with‐migration model. Mol Ecol Resour, 19, 1593–1609.
Hudson,R.R. and Kaplan,N.L. (1985) Statistical properties of the number of recombination events in the history of a sample of DNA sequences. Genetics, 111, 147–164.
Jennings,W.B. and Edwards,S.V. (2005) Speciational history of Australian Grass Finches (Poephila) inferred from thirty gene trees. Evolution, 59, 2033–2047.
Mirarab,S. et al. (2014) ASTRAL: genome-scale coalescent-based species tree estimation. Bioinformatics, 30, i541–i548.
Mirarab,S. et al. (2016) Evaluating summary methods for multilocus species tree estimation in the presence of incomplete lineage sorting. Syst Biol, 65, 366–380.
Rozas,J. et al. (2017) DnaSP 6: DNA Sequence Polymorphism Analysis of Large Data Sets. Molecular Biology and Evolution, 34, 3299–3302.
Yang,Z. (2015) The BPP program for species tree estimation and species delimitation. Current Zoology, 61, 854–865.

