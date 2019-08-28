# Adaptation
Analysis scripts for adaptation WGS data pre-processed by [breseq](http://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/)    
Deatherage, D.E., Barrick, J.E. (2014) Identification of mutations in laboratory-evolved microbes from next-generation sequencing data using breseq. Methods Mol. Biol. 1151: 165–188.    
PMID: [24838886](https://www.ncbi.nlm.nih.gov/pubmed/24838886)     

-----

## STEP 1: Copy and rename the required files from aerobio output directory to local directory.     

*Note: this step is specific for TvO lab server ```prince```. If you're not working on this server, you will need to assemble your .gd files another way. If you ARE working on this server, please navigate to ```/store/data/share/Adaptation``` before running any of these scripts.*     

Running aerobio/breseq produces a directory under ```ExpOut/[ExpID]/``` where ```[ExpID]``` is the aerobio experiment ID (e.g. ```170509_NS500751_0035_AHGTWJBGX2```).     
Under this directory, the ```Out``` subdirectory contains a separate directory for each sample. And under each sample-specific directory, there exists a machine-readable ```output.gd``` file. This is the file we need to run all downstream analysis. (N.b. The ```summary.html``` file is much more human friendly and readable for a first-pass view of the results). An example of a sample-specific sub-directory (named ```T4VanPopulation4```) is shown below:    
```
└── T4VanPopulation4
    ├── 01_sequence_conversion
    ├── 02_reference_alignment
    ├── 03_candidate_junctions
    ├── 04_candidate_junction_alignment
    ├── 05_alignment_correction
    ├── 06_bam
    ├── 07_error_calibration
    ├── 08_mutation_identification
    ├── data
    └── output
        ├── calibration
        ├── evidence
        ├── index.html
        ├── log.txt
        ├── marginal.html
        ├── output.done
        ├── output.gd
        └── summary.html
```
Since all output files are named "output.gd" it is also necessary to rename these copied files as their corresponding sample names. 
Use the "copy_gd.py" script to do this:    

```python copy_gd.py -e 190109_NS500751_0113_AHFCW3AFXY/ -o gd/testgd```    

(make sure to use your own aerobio ID, and a reasonable output directory- preferably under "gd")     

## STEP 2: Filter mutations based on frequency

This step is for removing low frequency mutations and background mutations (genotypic changes either already present in the parental strain, or mutations that are condition-independent - i.e. mutations that are observed in a control adaptation experiment).     
For this step, you need to define a number of things:     
* input directory that contains all files (use ```./``` which is the current directory, unless there are changes to the directory structure here)
* an comparison sheet (.csv) that specifies which gd files to use as control and which to use as experiment    
* a genome file (.gbk) - genomes we have used for TIGR4, 19F and D39 are under the ```gbk``` directory    
* output file name (.csv) for filtered results    
* lower cutoff (```L```) - mutations with frequency >```L``` in the CONTROL experiments will be removed    
* upper cutoff (```U```) - mutations with frequency <```U``` in the EXPERIMENTAL conditions will be removed 
* optional: if the samples in the Comparison Sheet are **clonal** use the ```-C``` flag to run in clonal mode
The comparison sheet needs to be comma-separated, with 2 columns titled ```File``` and ```Group```. The ```File``` column lists the filename and location of the gd files, and the ```Group``` column has values ```E``` or ```C``` depending on whether that gd file is experiment or control.      
    
Example comparison sheets are in the directory ```ComparisonSheets```     
Example use of the ```filter_gd.py``` script in population mode:    

```python filter_gd.py -i ./ -s ComparisonSheets/T4_VNC.csv -g gbk/NC_003028.gbk -o filtered/T4_VNC_10_50.csv -l 10 -u 50```    

This will generate a file ```filtered/T4_VNC_10_50.csv``` that has aggregated the population frequencies of each mutation (rows) in each of the experimental or control conditions (columns). The last column contains the locus tag, if the mutation coordinate is contained within a coding sequence.    
    
Example use in clonal mode:    
```python filter_gd_test.py -i ./ -s ComparisonSheets/Suyen_RIF3_clone.csv -g gbk/NC_003028.gbk -o filtered/Suyen_RIF3_clone.csv -l 10 -u 50 -C```    
**NOTE**: if you want to aggregate all results from a set of experiments WITHOUT filtering out any mutations, set L=100 and U=0 (this might take longer)    

As an example, I did the filtering for the 6 ABX experiments in T4    

```
python filter_gd.py -i ./ -s ComparisonSheets/T4_CIP.csv -g gbk/NC_003028.gbk -o filtered/T4_CIP_10_50.csv -l 10 -u 50    
python filter_gd.py -i ./ -s ComparisonSheets/T4_KAN.csv -g gbk/NC_003028.gbk -o filtered/T4_KAN_10_50.csv -l 10 -u 50    
python filter_gd.py -i ./ -s ComparisonSheets/T4_LVX.csv -g gbk/NC_003028.gbk -o filtered/T4_LVX_10_50.csv -l 10 -u 50    
python filter_gd.py -i ./ -s ComparisonSheets/T4_PEN.csv -g gbk/NC_003028.gbk -o filtered/T4_PEN_10_50.csv -l 10 -u 50    
python filter_gd.py -i ./ -s ComparisonSheets/T4_RIF.csv -g gbk/NC_003028.gbk -o filtered/T4_RIF_10_50.csv -l 10 -u 50    
```
    
#### Generating input table for muller plots    
If you wish to use the [muller](https://github.com/cdeitrick/muller_diagrams/) package for downstream analysis of your populations, you can generate a "muller-formatted" filtered table. This will have all necessary columns for the muller package to process the experiment.    
Note: When setting up your comparison sheet, the Experimental samples should be the different timepoints of **the same population**. And the timepoints can be named in any way, as long as there is a **single** integer in the sample name (e.g. "Bio_day_12" or "LVX36" etc). The columns corresponding to each timepoint will be renamed with the integer contained in the sample name. You can use the ```filter_gd.py``` script in "muller" mode by adding the optional argument ```-M```. This will generate a second ```.csv``` file that is marked ```[outputfilename]_Muller.csv```, which is appropriately formatted for downstream use with the muller package.     
Example usage:     
```python filter_gd.py -i ./ -s ComparisonSheets/Biofilm_LVX_forMuller.csv -g gbk/NC_003028.gbk -o filtered/Biofilm_LVX_forMuller_10_30_test.csv -l 10 -u 30 -M```    


## STEP 3: Agregate results from multiple experiments, generate plots. 

For this, you need to define an experiment sheet. This will be a .csv file with three columns: ```File```, ```Name```, and ```MergeBy```. ```File``` will be the full file name of the filtered csv file obtained from Step 2. ```Name``` can be whatever descriptive name you choose to name that experiment. ```MergeBy``` is the column name as it appears on the ```Annotation_3Strains_Clean.csv``` metadata/annotation table to merge the locus tags by. Use the following:     
* T4: "TIGR4.old"    
* 19F: "TAIWAN.19F.new"    
* D39: "D39.new"    

Use the ```expt_analysis.py``` script:     

```python expt_analysis.py -i ExperimentSheets/T4_ABX.csv -o Results/T4_ABX```    

This will generate a directory "T4_ABX" under "Results", and save a number of files in here.     
* ```Mutation_summary.csv```: This is a summary of the number of mutations and the types of mutations (INS/DEL/SNP/SUB; Coding/Noncoding), Total number of mutations, and the total number of unique locus tags with mutations (which is the number of adapted genes)
* ```Mutation_type.svg```: Bar plot of numbers of mutations of each type (SNP/INS/DEL/SUB) in each experimet    
* ```Mutation_coding.svg```: Bar plot of numbers of mutations in Coding vs Noncoding sequences in each experiment    
* ```Mutation_annotation_all.csv```: each adapted gene from each experiment is aggregated, and merged with its annotation  
* ```AG_TAG_mutations.svg```: Bar plot of the Functional Tags of adaptive mutations in each experiment    
* ```AG_CATEGORY_mutations.svg```: Bar plot of the Functional Categories of adaptive mutations in each experiment  
* ```AG_TAG_genes.svg```: Bar plot of the Functional Tags of adapted genes in each experiment    
* ```AG_CATEGORY_genes.svg```: Bar plot of the Functional Categories of adapted genes in each experiment    

