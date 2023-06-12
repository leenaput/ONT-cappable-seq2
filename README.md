# ONT-cappable-seq2

We describe here a general overview of the updated ONT-cappable-seq data analysis pipeline, which is now fully automated. For a more detailed overview of individual steps, see the old ONT-cappable-seq repository.

### **I. Data navigation**

The repository is organized in the following way:

ONT-cappable-seq2/
- workflow/
	- envs/
	- rules/
	- snakefile
- config/
	- config.yaml
- input/
	- fastq_data/
	- genome_data/
- peak_clustering.r

### **II. Setting up the configuration file**

The config.yaml file needs to be modified to your experimental settings. Adjust the paths to your input files and the different parameters used for annotation of the transcriptional boundaries. 

where:
- **Sample name**:  contains the name of your organism of interest (sampleID, *e.g LUZ19*)
- **Fasta file**: contains the path to your reference genome of interest (input/genome_data/sampleID.fasta, *e.g input/genome_data/LUZ19.fasta*)
- **Enriched fastq**: contains the path to the raw sequencing file of your enriched sample (input/fastq_data/sampleID_enriched.fastq, *e.g input/fastq_data/LUZ19_enriched.fastq*)
- **Control fastq**: contains the path to the raw sequencing file of your control sample (input/fastq_data/sampleID_control.fastq, *e.g input/fastq_data/LUZ19_control.fastq*)
- **ID**: additional identifier to classify your samples (cannot be empty), e.g timepoint, specific condition, replicate, etc.
- **Termseq alpha**: peak calling threshold value used by the termseq-peak algorithm
- **Cluster width**: peak positions within the specified distance are clustered. The position with the highest number of reads is taken as the representative of the cluster.
- **Minimum coverage**: absolute number of reads required at this position to be considered as candidate TSS/TTS position. This can be specified for the enriched and the control sample separately. 
- **Peak alignment error**: Positional difference (n) allowed between peak positions identified in the enriched and the control dataset used to calculate the enrichment ratio at a specific genomic position i, based on the read count per million mapped reads (RPM) at that peak position, as   calculated by: 

			enrichment ratio (i)= (RPM_(enriched sample) (i))/(RPM_(control sample) (i±n)) 
- **TSS threshold**: enrichment ratio value that needs to be surpassed to annotate 5’ peak position as a TSS.
- **TTS threshold**: minimum read reduction to annotate 3’ peak position as TTS, determined calculating the coverage drop across the putative TTS, averaged over a 20bp window
- **TSS/TTS sequence extraction**: selection of promoter and terminator region for which you want to extract the DNA sequence, defined by the number of up- and downstream nucleotides relative to the TSS and TTS.
