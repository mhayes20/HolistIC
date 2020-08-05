# HolistIC: Unambiguous prediction of multiple double minute chromosome architectures
```
Usage: python holistic.py <in.amplicons.bed> <in.gothic.csv> <float: q-value threshold> <autosomes only Y|N>
		in.amplicons.bed:      BED file containing coordinates of identified eccDNA segments (amplicons). Column 4 should uniquely identify the contig to which the amplicon belongs.
		in.gothic.csv:         CSV file containing GOTHiChicup output dumped from R data frame.
		q-value threshold:     Count interaction as significant if q-value is less than or equal this value. (Recommended: 0.05).
		autosomes?:            If Y, HolistIC will look for interactions on the whole genome including X and Y chromosomes.
```

## Input
### BED file

The BED file must be a sequence of records of the following format:
```
chr	start	end	ID
```

Assume that we have the following NGS-based predicted double minutes (i.e. where we know the breakpoints of the amplicons and their adjacencies):

<img width="594" alt="holistic github fig1 " src="https://user-images.githubusercontent.com/10326087/89368209-8c966780-d6a0-11ea-8da9-4eda862993e2.png">

The BED file will contain the following records:

```
chr8	28000000	28500000	DM1
chr19	40000000	40200000	DM1
chr5	3000000		3300000		DM1
chr4	2000000		2300000		DM2
chr11	6000000		6800000		DM2
chr18	20000000	20200000	DM2

```


They will have the following BED input file:

The BED file can be created using the output from programs that reconstruct circular eccDNA (like AmpliconArchitect or CouGaR).

### CSV file


## Output

Assume that we have the following predicted double minute chromosome:

```



```
The output has the following format:
