# miRTranscript

<pre>
miRTranscript is an open source python implementation made available under the GNU General Public License. It requires Python3 or above for smooth running. The software screens novel precursor sequences from transcriptome datasets<br> of plants. There is a huge inflow of data in terms of genomics and molecular biology. There are several softwares for screening miRNA from small RNA-seq data. However, there is a scarcity of working softwares for screening miRNA from transcriptome data. Moreover the ones available are mostly based on homology. Here, we present an ML-based system with Random Forest algorithm named miRTranscript for the computational prediction of miRNA from plant transcriptome datasets. 

If you're using miRTranscript please cite us:
Ayachit, G., Pandya,H., Das, J. (2020).miRTranscript: Pipeline for automated detection of novel miRNAs from plant transcriptomes using a mixed homology and machine learning approach.Manuscript submitted for publication. 

System requirements:
-Ubuntu (>=16.04)(The current version does not yet support other platforms or run via conda environment)

Dependencies required to run miRDetect:
-Blast executables
-ViennaRNA Package 
-Biopython
-scikit-learn


Building local server
git clone https://github.com/Garima268/miRTranscript.git
cd miRTranscript
python install.py
Download uniprot_sprot fasta from https://www.uniprot.org/downloads and place the file in DB folder

Once the dependencies are installed please enter the absolute paths of all the above in the config.py 
For example
  #Enter full path to blast executables
  path_blast = "/usr/bin/ncbi-blast-2.10.0+/bin"
  #Enter full path to ViennaRNA package directory
  path_vienna = "/usr/bin/ViennaRNA-2.4.14/"
  #Enter full path to Uniprot/nr database directory
  path_db = "/home/User/miRDetect/DB"
  #Enter name of Database fasta file
  p_name = "uniprot_sprot.fasta"
  

Before running make sure that the input fasta file headers do not have spaces in them. Please remove spaces in headers from file

Usage:
python miRTranscript.py [options] Assembled_transcriptome_fasta 
Options:
-p int Specify Number of threads to be used for Blast 

</pre>

