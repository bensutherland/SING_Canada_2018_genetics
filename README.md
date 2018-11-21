# SING_Canada_2018_genetics
Associated code for analysis of SING 2018 genetic data.    
More information about the Summer Internship for INdigenous Peoples in Genomics (SING) Canada program can be found here: http://indigenoussts.com/sing-canada/sing-canada-2018/

#### Requirements:    
`vcftools`    https://vcftools.github.io/index.html    
`adegenet`    https://cran.r-project.org/web/packages/adegenet/index.html    
`vcfR`    https://cran.r-project.org/web/packages/vcfR/index.html    
`hierfstat`     https://cran.r-project.org/web/packages/hierfstat/index.html    
`ape`     https://cran.r-project.org/web/packages/ape/index.html    
`jvarkit`    http://lindenb.github.io/jvarkit/MsaToVcf.html    

`java jdk 8.0`    http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
To use jvarkit, you will need java jdk 8.0. Check your java: 
Check java:
`/usr/bin/java -version`     

#### Inputs:
multiple sequence alignment file (msa)    

### 1. Prepare alignment file into a VCF
Put the .aln file into `02_raw_data`. In this example, the file is entitled `BaynesSouthern_Butter.aln`          

Move from the sequence alignment file into a vcf file (provide path to msa2vcf)    
`cat 02_raw_data/BaynesSouthern_Butter.aln | /Users/wayne/programs/jvarkit/dist/msa2vcf --allsites > 03_vcf/aln.vcf`
The inclusion of `--allsites` is important as it will give you readouts for all loci.   

### 2. Filter the VCF
How many SNPs were identified before filtering?
`grep -vE '^#' 03_vcf/aln.vcf | wc -l`
111 sites are viewed in the vcf. (note this will not work with the --allsites flag, as it will include all loci in the VCF.      

And how many samples?    
`grep '|' 02_raw_data/butterclam.aln  | awk '{ print $1}' - | sort -n | uniq  | wc -l`    
84 samples are present in the vcf after removing the other species.   

With this dataset, there are 84 individual samples. If we want to make sure that at least 2 individuals have the SNP, then require 4 copies (because this is a diploid VCF with a haploid dataset, each individual is represented by 2 copies).  

Use VCFtools to filter the VCF and only keep sites that are present at least four times (2 individuals). For our data, this would equate to a minor allele frequency filter (MAF) of 4 / (84 * 2) = 0.0238 (or MAF > 0.02).   
`vcftools --vcf 03_vcf/aln.vcf --recode --maf 0.024 --out 03_vcf/aln_maf`     
Keeps 4 of 111 sites.   

The SNPs with the highest MAF is at 400 bp, but remember that these are all in *tight physical linkage* so technically we should only be using a single marker per Fst evaluation. In this case, there are only five markers retained, and so this is probably fine to just keep together. If necessary, we could just retain the highest MAF polymorphism.       
In fact, if there was an option to do error correction to prevent false-positive polymorphisms on this data, the entire marker itself should be treated as a haplotype marker and this could be used as a single marker input for Fst evaluation.      
In any case, the results provided here most likely provide the same answer as that would be obtained from such an analysis.    

### 3. Generate stats on data (Fst)
Use the Rscript `01_scripts/vcf_to_hfstat.R` to turn the vcf into a genlight file, then assign population names to the samples based on sample names, plot a neighbour-joining tree, translate the data to hierfstat format and calculate Fst. This will also produce a graph of sample size per population.    
Outputs include:    
* `04_results/njt.pdf` (neighbour-joining tree)
* `04_results/sample_size_per_pop.pdf` (barplot of sample sizes)
* `04_results/all_data_wcfst.csv` (Weir Cockerham Fst)
* `04_results/fst_nboot_<#bootstraps>.csv` (WC Fst with bootstraps for significance)
