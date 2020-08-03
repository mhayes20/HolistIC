# HolistIC: Unambiguous prediction of multiple double minute chromosome architectures
Quick start:
Usage: python holistic.py <in.amplicons.bed> <in.gothic.csv> <float: q-value threshold> <autosomes only Y|N>
in.amplicons.bed:	BED file containing coordinates of identified eccDNA segments (amplicons). Column 4 should uniquely identify the contig to which the amplicon belongs.
in.gothic.csv:		CSV file containing GOTHiChicup output dumped from R data frame.
	q-value threshold:	Amplicon interaction is significant if q-value of bins' interaction is less than or equal this value. (Recommended: 0.05).
	autosomes?:	If Y, HolistIC will look for interactions on the whole genome including X and Y chromosomes.
 
