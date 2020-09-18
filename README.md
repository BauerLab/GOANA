### GOANA

GOANA is a tool for assessing mutation rates at target sites across the genome from next generation sequencing data.
 
As input it takes alignment files (in bam format) for a control and treated sample as well as information about the position of the target sites provided in bed format (chr, start, stop as a minimum). The target sites are compared across the two samples and the mutation rate calculated. Users can specify thresholds and filters to reduce the impact of low-frequency mutations or SNPs.
 
GOANA will generate a single report (.txt format) for all target sites.
 
A web-app implementation is also available from gt-scan.csiro.au/goana.

---
**Usage**  
```
GOANA.py [-h] [-s] [-q [Q]] [-up [UP]] [-down [DOWN]] [-mc [MC]] [-mr [MR]] [-o [O]] regions.bed control.bam treated.bam
```
---

__positional arguments:__  
&nbsp;&nbsp;&nbsp;&nbsp;__regions:__ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      List of regions (bed)  
&nbsp;&nbsp;&nbsp;&nbsp;__control:__ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      Control BAM file  
&nbsp;&nbsp;&nbsp;&nbsp;__treated:__ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;      Treated BAM file  
  
__optional arguments:__  
&nbsp;&nbsp;&nbsp;&nbsp;__-h, --help:__ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  	 	Show this help message and exit  
&nbsp;&nbsp;&nbsp;&nbsp;__-up [UP]:__   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;		Extend all regions upstream by X positions  
&nbsp;&nbsp;&nbsp;&nbsp;__-down [DOWN]:__ 	&nbsp;&nbsp;&nbsp;&nbsp;Extend all regions downstream by X positions  
&nbsp;&nbsp;&nbsp;&nbsp;__-mc [MC]:__  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    		Minimum read coverage over region (decimal between 0 and 1)  
&nbsp;&nbsp;&nbsp;&nbsp;__-mr [MR]:__  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   		Minimum relative read frequency to classify as significant mutation (percentage between 0 and 100)  
&nbsp;&nbsp;&nbsp;&nbsp;__-o [O]:__
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Optional output file   
&nbsp;&nbsp;&nbsp;&nbsp;__-s:__
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Include this flag to exclude ALL SNPs

&nbsp;&nbsp;&nbsp;&nbsp;__-q [Q]:__
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Only include bases with a quality (Phred score) equal to or above this threshold (INT)

  
__Requires:__  
&nbsp;&nbsp;&nbsp;&nbsp;Samtools + Pysam  
&nbsp;&nbsp;&nbsp;&nbsp;Numpy  
&nbsp;&nbsp;&nbsp;&nbsp;Python3.7.2
