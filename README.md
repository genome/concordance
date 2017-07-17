# concordance
CWL and py script for tool to report on concordance between tumor and normal samples

Command: python newConcordance.py bam1 bam2 ref_fasta snps output --output_geno geno

(output_geno is optional, will create a list of all genotypes determined by the tool)

**this tool's CWL can be found in cancer-genomics-workflow and can be run in Docker container mgibio/docker-concordance**
