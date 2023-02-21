# BarBIQ Pipeline (v1.2.0)

The BarBIQ pipeline is for processing the sequencing data (please note that this code only accept the phred+33 encoding system for the sequencing quality scores) obtained by BarBIQ method to identify Bar sequences (16S rRNA sequences) and cOTUs (cell-based taxa) and to quantify the identified cell numbers of each cOTU. 

This is version 1.2.0\
Author: Jianshi Jin\
Supervisor: Katsuyuki Shiroguchi<br/>


## What's new in this version compared to version 1.1.0

(1) In BarBIQ_Step_2.pl, the lowest-level taxon of each cOTU predicted by the RDP classifier with a classification bootstrap confidence equal or higher than a given threshold (between 0 and 1, which can be set up in the BarBIQ_example_parameter_2.txt file before running BarBIQ_Step_2.pl) is printed as the last column of the "Absolute-concentration" sheet in BarBIQ_Quantification.xlsx; for details about the classification bootstrap confidence in RDP classifer, please refer to https://github.com/rdpstaff/classifier and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965039/.<br/> 

(2) An optional process is added. For changing the threshold of classification bootstrap confidence to print the lowest-level taxon (see above (1)) after the processing of BarBIQ_Step_2.pl, the code BarBIQ_final_bac_con_nomalization_xlsx.pl can be run independently with a given threshold (between 0 and 1) (an instruction to use this code is newly provided below).<br/>

(3) In BarBIQ_Step_2.pl, a new file is output and named as BarBIQ_Bar_sequences.fasta. Bar sequences are printed in fasta format with taxonomic labels (predicted by the RDP classifier) and cOTU ID which may be used for other downstream bioinformatic analyses.<br/>

(4) A size-reduced demo data is provided for the quick testing of this pipeline; the instruction to use this demo data is provided below.<br/>

(5) In BarBIQ_Step_1.pl, the processed data are rearranged for saving memory.<br/>

## Contents
(1) BarBIQ_code: All codes of this pipeline<br/>

(2) Demo: Example data which can be used for demonstration.<br/>
* There are two types of data. The "Data_and_inputfiles_for_whole_processing" data can be used to successively run both BarBIQ_Step_1.pl and BarBIQ_Step_2.pl, but it requires 6 + 0.5 hours of running time; the "Data_and_inputfiles_for_quick_testing" data can be used to quickly test BarBIQ_Step_1.pl in 1 min and BarBIQ_Step_2.pl in 5 min, respectively (we note that this data cannot be used to successively run both BarBIQ_Step_1.pl and BarBIQ_Step_2.pl). <br/>

(3)Expected_output_files: Expected output after this pipeline runs correctly<br/>

## Codes
The codes of this pipeline is written in Perl (v5.22.1) or in R (version 3.6.3)(only one code). <br/>

## System requirements
(1) Required operation system: Linux 64<br/>

(2) Required software: Perl(tested v5.22.1); R(tested versions 3.5.1, 4.0.2, and 4.1.1); nucleotide-sequence-clusterizer (Version 0.0.7)<br/>

(3) Required perl modules: IPC::System::Simple; File::Path; File::Copy; Parallel::ForkManager; Bio::SeqIO; Bio::Seq; Text::Levenshtein::XS; Text::WagnerFischer; List::Util; Statistics::Basic; Excel::Writer::XLSX; Math::Round; Bundle::BioPerl<br/>
* Perl modules can be installed by the code BarBIQ_check_and_install_modules.pl (see below).<br/>

(4) Required R packages: plotrix<br/>
* R packages can be installed by the code BarBIQ_final_fitting_OD.r during the processing of BarBIQ_Step_2.pl.<br/>


## Download

Clone/download this repository to your local computer<br/>
Click "Download ZIP" button<br/>
Or<br/>
Type in the Terminal or equivalent:<br/> 
```
$ git clone https://github.com/Shiroguchi-Lab/BarBIQ_Pipeline_V1_2_0.git
```

## Install required software if need
(1) Install Perl(https://www.perl.org)<br/>
(2) Install R(https://www.r-project.org)<br/>
(3) Install nucleotide-sequence-clusterizer (http://kirill-kryukov.com/study/tools/nucleotide-sequence-clusterizer/)<br/>
(4) Install RDPTools (https://github.com/rdpstaff/RDPTools)<br/>

## Install required perl modules
Type in the Terminal or equivalent with adjustment of the path, "path-to", accordingly:
```
$ cd path-to/BarBIQ_code
$ chmod a+x *
$ ./BarBIQ_check_and_install_modules.pl
```

## How to use the codes for real data or the demo data "Data_and_inputfiles_for_whole_processing"
### Step 1: For each sequencing run
(A) Modify the "Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_inputfile_1.txt" and "Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_parameter_1.txt" for your own data (make sure to adjust the paths, "path-to", accordingly, even for the demonstration)<br/>
* Details are explained in the file itself or in Instructions_prepare_input_files.txt<br/>

(B) Type the following command with adjustment of the path, "path-to", accordingly:<br/>
  ```
$ path-to/BarBIQ_code/BarBIQ_Step_1.pl --in path-to/Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_inputfile_1.txt --p path-to/Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_parameter_1.txt
  ```

### Step 2: For all samples
(A) Modify the "Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_inputfile_2.txt" and "Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_parameter_2.txt" for your own data (make sure to adjust the paths, "path-to", accordingly, even for the demonstration)<br/>
* Details are explained in the file itself or in Instructions_prepare_input_files.txt<br/>

(B) Type the following command in the Terminal or equivalent with adjustment of the path, "path-to", accordingly:<br/>
  ```
$ path-to/BarBIQ_code/BarBIQ_Step_2.pl --in path-to/Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_inputfile_2.txt --p path-to/Demo/Data_and_inputfiles_for_whole_processing/BarBIQ_example_parameter_2.txt
  ```

## How to quickly test the codes using the demo data "Data_and_inputfiles_for_quick_testing"
### For Step 1
(A) Adjust the paths, "path-to", accordingly in the "Demo/Data_and_inputfiles_for_quick_testing/For_Step1/BarBIQ_example_parameter_1_S.txt"<br/>

(B) Type the following command with adjustment of the path, "path-to", accordingly:<br/>
  ```
$ path-to/BarBIQ_code/BarBIQ_Step_1.pl --in path-to/Demo/Data_and_inputfiles_for_quick_testing/For_Step1/BarBIQ_example_inputfile_1_S.txt --p path-to/Demo/Data_and_inputfiles_for_quick_testing/For_Step1/BarBIQ_example_parameter_1_S.txt
  ```

### For Step 2
(A) Adjust the paths, "path-to", accordingly in the "Demo/Data_and_inputfiles_for_quick_testing/For_Step2/BarBIQ_example_inputfile_2_S.txt" and "Demo/Data_and_inputfiles_for_quick_testing/For_Step2/BarBIQ_example_parameter_2_S.txt"<br/>

(B) Type the following command in the Terminal or equivalent with adjustment of the path, "path-to", accordingly:<br/>
  ```
$ path-to/BarBIQ_code/BarBIQ_Step_2.pl --in path-to/Demo/Data_and_inputfiles_for_quick_testing/For_Step2/BarBIQ_example_inputfile_2_S.txt --p path-to/Demo/Data_and_inputfiles_for_quick_testing/For_Step2/BarBIQ_example_parameter_2_S.txt
  ```

## Optional analysis
To print the lowest-level taxon of each cOTU predicted by the RDP classifier with a classification bootstrap confidence equal or higher than a given threshold (between 0 and 1; using "--cutoff 0.X" to adjust) after the processing of BarBIQ_Step_2.pl, by using the code BarBIQ_final_bac_con_nomalization_xlsx.pl

  ```
$ cd path-to/Output_Step2
$ path-to/BarBIQ_code/BarBIQ_final_bac_con_nomalization_xlsx.pl --count *_annotation_RDP_classifier.txt --total Total_concentration.txt --out BarBIQ_Quantification.xlsx --taxa Classifier --cutoff 0.8
  ```

## Key output results
### For Step 1
#### (1) Filename end with "_clean_sub1_sub2_shift_link_derr_indels_chimeras_similar_statistic"
* In the output folder of each index of the Step 1 contains the number of droplets (2nd column) for each identified Bar sequence (3rd column); in the 4th column, "LINK" and "I1R2" indicates I1 and R2 are overlapped or not; the numbers in the 5th and 6th columns indicate the numbers of bases obtained from I1 and R2; the numbers in the 1st column are serial number.<br/> 
* For the detailed information, please refer to the "Step 10" in Supplementary Note 2 of [High-throughput identification and quantification of single bacterial cells in the microbiota](https://www.nature.com/articles/s41467-022-28426-1).<br/>

#### (2) Filename end with "_clean_sub1_0.1_sub2_shift_link"
* In the output folder of each index of the Step 1 contains the "SCluster" ID, RepSeqs for I1 reads, RepSeqs for R2 reads, "LINK" or "I1R2" indicating I1 and R2 are overlapped or not, and linked RepSeqs, which is used for identifying multiple Bar sequences from the same bacterium. <br/>
* For the detailed information, please refer to the "Step 15" in Supplementary Note 2 of [High-throughput identification and quantification of single bacterial cells in the microbiota](https://www.nature.com/articles/s41467-022-28426-1).<br/>

#### (3) Filename with "_reads_per_barcode_aveReads_" in the middle
* In the output folder of the Step 1, there are one file per one index; each file contains the detected Barcodes and the number of detected reads for each barcode; the number in the filename after "_aveReads_" indicates the average number of reads per barcode, and the number in the end of the filename indicates the detected number of barcodes.<br/>

### For Step 2
#### (1) "BarBIQ_Quantification.xlsx" 
* In the output folder of the Step 2 contains the raw counted cell numbers, absolute cell numbers per unit weight or volume, and proportional cell abundance of the cOTUs for each sample with the taxonomies predicted by RDP classifier. <br/>
* For the detailed information, please refer to the "ReadMe" sheet in "BarBIQ_Quantification.xlsx". <br/>

#### (2) "BarBIQ_Bar_sequences.xlsx" 
* In the output folder of the Step 2 contains the Bar sequences for all samples with the cOTU information. <br/>
* For the detailed information, please refer to the "ReadMe" sheet in "BarBIQ_Bar_sequences.xlsx". <br/>

#### (3) "BarBIQ_Bar_sequences.fasta"
* In the output folder of the Step 2 contains Bar sequences in fasta format with taxonomic labels (predicted by the RDP classifier) and cOTU ID. <br/>
* This file may be used for other downstream bioinformatic analyses such as phylogenetic tree analysis. <br/>

## Documentation
All the algorithms of this pipeline is described in the paper [High-throughput identification and quantification of single bacterial cells in the microbiota](https://www.nature.com/articles/s41467-022-28426-1).<br/>

## How to cite the BarBIQ pipeline
### Please cite the following two papers:
Jin, J., Yamamoto, R., Takeuchi, T., Cui, G., Miyauchi,E., Hojo, N., Ikuta, K, Ohno, H. & Shiroguchi, K. High-throughput identification and quantification of single bacterial cells in the microbiota. Nat. commun., 13, 863, (2022). <br/>

Ogawa, T., Kryukov, K., Imanishi, T. & Shiroguchi, K. The efficacy and further functional advantages of random-base molecular barcodes for absolute and digital quantification of nucleic acid molecules. Sci. Rep. 7, 13576 (2017).<br/>

### If you use the taxonomies for the cOTUs generated by this pipeline, please cite this paper as well:
Wang, Q, Garrity, G. M., Tiedje, J. M. & Cole, J. R. Naïve Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Appl. Environ. Microbiol. 73(16), 5261-5267 (2007).<br/>

## Other versions
[v1.0.0](https://github.com/Shiroguchi-Lab/BarBIQ)<br/>
[v1.1.0](https://github.com/Shiroguchi-Lab/BarBIQ_Pipeline_V1_1_0)<br/>
