# Leapfrog - the tool

Developer: Mark Fiers [https://github.com/mfiers/leapfrog](https://github.com/mfiers/leapfrog)

Data: Paired genomic sequence data

## Outline
Leapfrog attempts to take paired data and determine transposable element locations when compared to a reference sequence. Leapfrog takes paired data and maps to a repeat sequence database, alignments are filtered based on whether a read maps while its mate does not.  
The unmapped mate sequence is then extracted and mapped a reference sequence of particular species.  
The pileup of reads is then profiled for clusters which provide evidence of a repeat element. This is output a gff file.

## The current steps of leapfrog

### lf_danglers

lf_danglers uses the Bowtie2 aligner to map paired reads to a repeat sequence library.   
Information about mapping of either read in a pair is checked using samtools bit scores.   
Currently no attempt to identify junction spanning reads, just 'danglers' which are unmapped reads with a mate that maps to the repeats.   
Information about the mapped mate is captured in the fastq of unmapped reads in the naming, this is carried through to the next step.

### Mapping to the reference sequence

Currently the user must do this step themselves.  
The expected output for the next step is from Bowtie2 as aligner again. This time reads are mapped to the genome reference as single end reads.  

### lf_regionify

lf_regionify takes the bamfile output after aligning to the reference sequence.  
Potentially some flags that are unique to Bowtie2 use but this needs checking.  
Clusters 'peaks' are attempted to eb identified by coverage plot changes.   
A minimal coverage is input by the user with a defualt value otherwise.  
Peaks are distinguished by change differences which have a default but can also be set by the user.  
Estimates of average coverage are estimated by peak size and read length & count.  
An attempt at determining strand is also done, though could be improved and broadened. An estimation of a unique/nonuniqueness to the reads or insert are estimated.

### lf_finddiff 

lf_finddiff takes in a series of gff files and attempts to compare them. This is experimental according to the readme and has only had a brief development attempt.  
Currently compares entries in gff files based on location and family.  
Unsure how deals with slight difference in location estimation....

### lf_diffsum

lf_diffsum, again believe a draft was written but not thoroughly tested.  
Possibly trying to estimate overlap of inserts and generate a venn diagram?

## Possible improvements

1. Aligner used, try baw-MEM? Or some others and get some alignment stats for comparison
2. Generate simulated data for comparison and testing any modifications (have a simple pipeline and tool for this which can be developed easily)
3. Attempt to extract junction spanning reads
4. Check and modify filtering of danglers and regionify steps:
	- improve the extraction
	- to extract some reports at the various stages
5. Pipeline the steps better and increase options for extraction and filtering
6. Check the regioinfy method, can this be improved?
7. Check the diff method
	- what is the variation in peak detection?
	- how can we compare peaks, strand, clustering etc
8. Develop some reporting on final outputs.




