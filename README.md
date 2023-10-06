**ExRec**

ExRec (Exclusion of Recombined DNA) is a Python pipeline that implements the four-gamete test to filter
out recombined DNA sites from up to thousands of DNA sequence loci. The pipeline consists of five 
standalone applications: the first two convert folders of NEXUS or PHYLIP files into the standard input file 
for the main program that conducts the four-gamete filtering procedures. The pipeline outputs 
recombination-filtered data in concatenated NEXUS and PHYLIP formats and a tab-delimited table 
containing descriptive statistics for all loci and the results. This software also allows the user to output 
the longest non-recombined sequence blocks from loci (current best practice) or randomly select non-
recombined blocks from loci (a newer approach). Two other applications in the package convert the 
recombination-filtered data into single-locus NEXUS or PHYLIP files. The ExRec package can thus facilitate 
species delimitation, species tree, and historical demography studies by providing loci that better meet 
the no-recombination assumption in coalescent-based analyses.

**Requirements and Installation**

The ExRec pipeline runs under Python 3 without any dependencies. You can install ExRec by 
downloading the five applications from the ExRec Github page. Each script is run using the 
command line.

**Documentation**

A user’s manual is also included in the ExRec package. This document contains the
instructions for using each of the five applications in ExRec as well as a tutorial that
uses two sample data sets (also included in the package). Instructions on how to
execute each application can also be found using the &quot;help&quot; command. On the
command line, just cd to the folder containing the application you want to run and then
type “help” on the command line (as shown below) followed by enter:

&gt;python3 Nexcombine.py help

&gt;python3 Phycombine.py help

&gt;python3 FGT.py help

&gt;python3 Nexsplit.py help

&gt;python3 Physplit.py help

License and Warranty

Please see the file &quot;LICENSE&quot; for details.
