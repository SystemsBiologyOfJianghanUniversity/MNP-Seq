This MNP-Seq software is used for the genotyping of amplicon sequecing datasets.
This software could be tested in linux envroments. The input could be in bam or 
fastq format. Both the Bowtie2 and BWA  sortware are required for the analysis, 
and you shoud install both of them before run MNP-seq.



Two steps should be conducted for the analysis, the first is analysed each pair of
reads based on the alingment with reference genome, this is completed by
genotyping4PE_step1.pl. And the final genotyping and comparision were finished
with genotyping4PE_step2.pl. If you have mutiple processors, it is better to
set up the parameter "-split" to use more processors to finish the jobs.
