
# HolistIC: Unambiguous prediction of multiple double minute chromosome architectures
```
Usage: python holistic.py <in.amplicons.bed> <in.gothic.csv> <float: q-value threshold> <autosomes only Y|N>
		in.amplicons.bed:      BED file containing coordinates of identified eccDNA segments (amplicons). Column 4 should uniquely identify the contig to which the amplicon belongs.
		in.gothic.csv:         CSV file containing GOTHiChicup output dumped from R data frame.
		q-value threshold:     Count interaction as significant if q-value is less than or equal this value. (Recommended: 0.05).
		autosomes?:            If Y, HolistIC will look for interactions on the whole genome including X and Y chromosomes.
```
## Prequisites
Tested with Python 2.7 but Python 3+ should work too</br>
Python Networkx: https://networkx.github.io/documentation/stable/install.html </br>
matplotlib: https://matplotlib.org/users/installing.html#installation-guide </br>
rpy2: https://pypi.org/project/rpy2/

## Input
The BED file are WGS-based coordinates of contiguous regions of copy number amplification, organized by contig. The CSV file is the dumped output from the GOTHiChicup function of the GOTHiC package. The input to HolistIC is data assumed to come from the same individual.

### BED file

The BED file must be a sequence of records of the following format:
```
chr	start	end	ID
```

Assume that we have the following NGS-based predicted double minutes (i.e. where we know the breakpoints of the amplicons and their adjacencies):

<img width="594" alt="holistic github fig1 " src="https://user-images.githubusercontent.com/10326087/89368209-8c966780-d6a0-11ea-8da9-4eda862993e2.png">

The BED file will contain the following records:

```
chr5	3000000		3300000		DM1
chr8	28000000	28500000	DM1
chr19	40000000	40200000	DM1
chr4	2000000		2300000		DM2
chr11	6000000		6800000		DM2
chr18	20000000	20200000	DM2
```

The BED file can be created using the output from programs that reconstruct circular eccDNA architectures (like AmpliconArchitect or CouGaR).

WARNING: Within each unique ID in the 4th column, the records must first be sorted in lexicographic order by chromosome, and then in ascending numeric order by the interval start position. This will not preserve the adjacencies provided in the BED file, but this does not adversely affect HolistIC's ability to segregate amplicons by contig. The input BED file can be sorted with the following command:

```
sort -k4,4 -k1,1 -k2,2 in.bed > out.bed
```

### CSV file
We assume that significant chromatin interactions are determined with GOTHiC. The CSV file is the dumped output from an R dataframe generated by the GOTHiChicup function:

https://master.bioconductor.org/packages/release/bioc/manuals/GOTHiC/man/GOTHiC.pdf

## Output

Assume that we have the following double minutes:

<img width="592" alt="Screen Shot 2020-08-04 at 11 41 58 PM" src="https://user-images.githubusercontent.com/10326087/89373196-ab4f2b00-d6ad-11ea-9754-cd557d9bc3c0.png">

Due to ambiguities inherent in WGS-only predictions, the above double minutes could have resulted in the following input BED file, where all amplicons are erroneously assumed to be part of the same double minute (remember that the records are lexicographically sorted by chromosome and then sorted ascending by start position of interval):

```
chr1	8000000		8900000		DM1
chr3	7000000		7400000		DM1
chr5	3000000		3300000		DM1
chr6	17000000	17800000	DM1
chr18	12000000	13000000	DM1
```
This is due to ambiguity induced by the red amplicon since the WGS data shows a red-blue-yellow-red-orange-green-red cycle.

If the CSV file of Hi-C interactions is provided, HolistIC will attempt to segregate each amplicon (i.e. row in BED file) into its true double minute. In the ideal case, HolistIC would have the following output for the above input BED file:
```
chr1	8000000		8900000		DM1
chr3	7000000		7400000		DM1
chr5	3000000		3300000		DM1
chr6	17000000	17800000	DM2
chr18	12000000	13000000	DM2
chr5	3000000		3300000		DM2
```

Furthermore, HolistIC will create PNG image files for each identified maximal clique, which ostensibly refer to separate double minutes. 
