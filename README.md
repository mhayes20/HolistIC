# HolistIC: Unambiguous prediction of multiple double minute chromosome architectures
```
Usage: python holistic.py <in.amplicons.bed> <in.gothic.csv> <float: q-value threshold> <autosomes only Y|N>
		in.amplicons.bed:      BED file containing coordinates of identified eccDNA segments (amplicons). Column 4 should uniquely identify the contig to which the amplicon belongs.
		in.gothic.csv:         CSV file containing GOTHiChicup output dumped from R data frame.
		q-value threshold:     Count interaction as significant if q-value is less than or equal this value. (Recommended: 0.05).
		autosomes?:            If Y, HolistIC will look for interactions on the whole genome including X and Y chromosomes.
```

## Input
### BED file format
The BED file must be a sequence of the records of the following format:
```
chr	start	end	ID
```

For example, for the following double minutes predicted in the same dataset:



They will have the following BED input file:


The BED file can be created using the output from programs that reconstruct circular eccDNA (like AmpliconArchitect or CouGaR).
