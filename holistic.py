import sys
import os.path
from os import path
import networkx as nx
from networkx import find_cliques
import matplotlib.pyplot as plt
from math import floor
import csv
import re
import networkx
from rpy2 import robjects
from rpy2.robjects.vectors import StrVector

chrs = ['chr1', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chrX', 'chrY']
#G = networkx.Graph() #Create the graph
firstChr = "chr1"
lastChr = "chrY"

g = nx.Graph()
pos=nx.spring_layout(g)
#G.add_nodes_from(["1 20 30", "50 2 20"])

dmCounter = 1
DMF_dict = {}

if(len(sys.argv) != 5):
	print ("Usage: python holistic.py <in.amplicons.bed> <in.gothic.csv> <float: q-value threshold> <autosomes only Y|N>")
	print("\t\tin.amplicons.bed:      BED file containing coordinates of identified eccDNA segments (amplicons). Column 4 should uniquely identify the contig to which the amplicon belongs.")
	print("\t\tin.gothic.csv:         CSV file containing GOTHiChicup output dumped from R data frame.")
	print("\t\tq-value threshold:     Count interaction as significant if q-value is less than or equal this value. Recommended: 0.05).")
	print("\t\tautosomes?:            If Y, HolistIC will look for interactions on the whole genome including X and Y chromosomes.") 
	sys.exit(0)


myfile = sys.argv[1]
 

myGothFile = sys.argv[2]
qthresh = float(sys.argv[3])
auto = sys.argv[4]

if (auto == "Y"):
	chrs = ['chr1', 'chr10','chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr2', 'chr20', 'chr21', 'chr22', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9']
	firstChr = "chr1"
	lastChr = "chr9"

dmfinder_data = []
goth_data = []
rb = 0
lb = 0
window_size = 200000



def importGothic():
	#This function will (somehow) import the data from Gothic

	#Will try to import into rpy2 dataframe 
	read_delim = robjects.r('read.delim')
	global goth_data
	global firstChr
	global lastChr

	goth_data = read_delim(myGothFile, header=False, sep=" ", stringsAsFactors = False)
	colnames = robjects.r('colnames')
	as_integer = robjects.r('as.integer')
	as_numeric = robjects.r('as.numeric')
	match = robjects.r.match

	
	goth_data.colnames = robjects.StrVector(['chrA','chrAstart','chrB','chrBstart','regArel','regBrel','prob','expected','readCount','pvalue','qvalue','logObsExp'])	
	my_col = match('chrAstart', goth_data.colnames)[0]
	goth_data[my_col - 1] = as_integer(goth_data[my_col - 1])
	my_col = match('chrBstart', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_integer(goth_data[my_col - 1])
	#Here: convert position columns to integer
	 
	my_col = match('regArel', goth_data.colnames)[0]
	goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])
	
	my_col = match('regBrel', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])

	my_col = match('prob', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])
	
	my_col = match('expected',goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])

	my_col = match('readCount', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_integer(goth_data[my_col - 1])
	
	my_col = match('pvalue', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])

	my_col = match('qvalue', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])

	my_col = match('logObsExp', goth_data.colnames)[0]
        goth_data[my_col - 1] = as_numeric(goth_data[my_col - 1])

	#print ("Length of Gothic data " + str(len(goth_data[:][0])))
	firstChr = goth_data[0][0]
	lastChr = goth_data[0][len(goth_data[:][0])-1]
	#print("NEXT CHR?: " + str(goth_data[1][8000]))
	#lastChr = goth_data[0][len(goth_data)-1]
	#print("LAST CHR AFTER GOTH IMPORT: " + str(lastChr))
	#print ("Will implement soon\n")

def countDM_DMFinder():
	#This function reads in the DMFinder output file and creates a dictionary where the keys are DM[ID integer] and the values are lists of DM amplicons and their coordinates

	colnames = robjects.r('colnames')
	#c = robjects.r('c')
	read_delim = robjects.r('read.delim')
	as_integer = robjects.r('as.integer')
	match = robjects.r.match
	global dmfinder_data
	dmfinder_data = read_delim(myfile, header=False, sep="\t", stringsAsFactors = False)
	#print (dmfinder_data)	
	

	dmfinder_data.colnames = robjects.StrVector(['chr','start','end','index', 'type'])
	#print ("Number of rows in data is " + str(len(dmfinder_data)))	
	my_col = match('start', dmfinder_data.colnames)[0]
	#print('Type of start before as.integer: %s' % dmfinder_data[my_col - 1].rclass[0])
	dmfinder_data[my_col - 1] = as_integer(dmfinder_data[my_col - 1])
	#print('Type of start after as.integer: %s' % dmfinder_data[my_col - 1].rclass[0])

	my_col = match('end', dmfinder_data.colnames)[0]
        #print('Type of end before as.integer: %s' % dmfinder_data[my_col - 1].rclass[0])
        dmfinder_data[my_col - 1] = as_integer(dmfinder_data[my_col - 1])
        #print('Type of end after as.integer: %s' % dmfinder_data[my_col - 1].rclass[0])
	
	
def bSearchFirstBin(position, start, end):
	
	global goth_data
	global window_size
	#print("position: " + str(position))
 	#print ("start: " + str(start))
        #print ("end: " + str(end))
	if(start >= end):
		#print "ERROR in bSearchFirstBin: value not found\n"
		return -1

	mid = int((end + start) // 2)

	target = position - (position % window_size)
	#print ("target: " + str(target))
	if(goth_data[1][mid] == target and goth_data[1][mid-1] != target):
		return mid
	elif(goth_data[1][mid] < target):
		return bSearchFirstBin(position, mid+1, end)
	else:
		return bSearchFirstBin(position, start, mid)
	
def bSearchChr(chr1, chr2, start, end):
	#Function performs a binary search, but the search key is the boundary of the two chromosomes

	global goth_data

	if(start >= end):
		#print "ERROR in bSearchChr: " + chr1 + " and " + chr2 + " boundary not found\n"
		return -1	#ERROR: boundary not found

	mid = int((end + start) // 2)
	#print ("Midpoint of binary search is: " + str(mid))
	#print ("chr1: " + chr1)
	#print ("chr2: " + chr2)
	#print ("start: " + str(start))
	#print ("end: " + str(end))
	#print("got_data[0][mid]: " + str(goth_data[0][mid]) + "\n")

	if(goth_data[0][mid] == chr1 and goth_data[0][mid+1] == chr2):
	#	print("A")
		return mid
	elif(goth_data[0][mid] <= chr1):
	#	print("B")
		return bSearchChr(chr1,chr2,mid+1,end)
	else:
	#	print("C")
		return bSearchChr(chr1,chr2,start,mid)
	
def bSearch2Chr(chr1, start, end): #DON'T THINK I"LL USE THIS
	#This binary search searches for the SECOND chromosome in an interaction within Gothic data. I.e. it searches within the chr1 parameter records and not throughout the whole data itself	
		
	global goth_data
	
	if(start >= end):
		return -1

	mid = floor((end - start)/2)

	if mid == 0:	#Prevent array out of bounds exception
		if goth_data[mid][2] == chr1:
			return mid
	else:		
		if goth_data[mid][2] == chr1 and goth_data[mid-1][2] < chr1 or goth_data[mid][2] == chr1 and goth_data[mid-1][0] != goth_data[mid][0]:	#Get the location of first record
			return mid
		elif goth_data[mid][2] >= chr1:
			return bSearch2Chr(chr1, start, mid)
		else:
			return bSearch2Chr(chr1, mid+1, end)
	
def chromBinSearch2(chromosome):
        #Search gothic data for the start and end position of the given chromosome
        global chrs
        global firstChr
        global lastChr
        global goth_data
        
        

        pos = chrs.index(chromosome)    #Get position in master chr list of the desired chromosome

	#print ("position of chr: " + chromosome + " is " + str(pos))
        #print ("Length of gothic data is: " + str(len(goth_data[0][:])))
	#print ("firstChr: " + firstChr)
        if chromosome == firstChr:      #If first chromosome in the master list, then it is the first chromosome in the records so zero is its index
                lb = 0 #left boundary is zero
                rb = bSearchChr(chrs[pos], chrs[pos+1], 0, len(goth_data[0][:])-1)                                 #Right boundary has to be searched
                if rb == -1:
                        print("1")
                        exit("No boundary found in GOTHiC data between chromosome " + chrs[pos] + " and " + chrs[pos+1] + "\n")

        elif chromosome == lastChr:     #Same as above but for the last record
                rb = len(goth_data[0][:]) - 1
                lb = bSearchChr(chrs[pos-1], chrs[pos], 0, len(goth_data[0][:])-1)
                if lb == -1:
                        print("2")
                        exit("No boundary found in GOTHiC data between chromosome " + chrs[pos-1] + " and " + chrs[pos] + "\n")
        else:                           #chromosome is somewhere in the middle so we have to search for its boundaries using two binary searches
                #Find first boundary    

                lb = bSearchChr(chrs[pos-1], chrs[pos], 0, len(goth_data[0][:])-1)
                if lb == -1:
                        print("3")
                        exit("ELSE STATEMENT: No boundary found in GOTHiC data between chromosome " + chrs[pos-1] + " and " + chrs[pos] + "\n")

                rb = bSearchChr(chrs[pos], chrs[pos+1], 0, len(goth_data[0][:])-1)
                if rb == -1:
                        print("4")
                        exit("ELSE STATEMENT: No boundary found in GOTHiC data between chromosome " + chrs[pos] + " and " + chrs[pos+1] + "\n")

        lb = lb + 2
        rb = rb + 2
        
        return (lb, rb)


def chromBinSearch(chromosome):
	#Search gothic data for the start and end position of the given chromosome
	global chrs
	global firstChr
 	global lastChr
	global goth_data
	global lb
	global rb	

	pos = chrs.index(chromosome)	#Get position in master chr list of the desired chromosome
	
	#print ("Length of gothic data is: " + str(len(goth_data[0][:])))
	#print ("position of chr: " + str(pos))

	#print ("lastChr: " + str(lastChr))
	
	#print ("position of chr: " + chromosome + " is " + str(pos) + " and firstChr is " + firstChr)
	if chromosome == firstChr:	#If first chromosome in the master list, then it is the first chromosome in the records so zero is its index
		lb = 0 #left boundary is zero
	#	print ("position of chr: " + str(pos))
                #print ("position of chr: " + str(len(chrs)))
		rb = bSearchChr(chrs[pos], chrs[pos+1], 0, len(goth_data[0][:])-1) 				#Right boundary has to be searched
		if rb == -1:
			print("1")
			exit("No boundary found in GOTHiC data between chromosome " + chrs[pos] + " and " + chrs[pos+1] + "\n")		
	
	elif chromosome == lastChr:	#Same as above but for the last record
	#	print ("position of chr: " + str(pos))
               # print ("position of chr: " + str(len(chrs)))
		rb = len(goth_data[0][:]) - 1
		lb = bSearchChr(chrs[pos-1], chrs[pos], 0, len(goth_data[0][:])-1)
		if lb == -1:
			print("2")
			exit("No boundary found in GOTHiC data between chromosome " + chrs[pos-1] + " and " + chrs[pos] + "\n")
	else:				#chromosome is somewhere in the middle so we have to search for its boundaries using two binary searches
		#Find first boundary	
			
	#	print ("position of chr: " + str(pos))
                #print ("position of chr: " + str(len(chrs)))
		lb = bSearchChr(chrs[pos-1], chrs[pos], 0, len(goth_data[0][:])-1)
		if lb == -1:
			print("3")
                        exit("ELSE STATEMENT: No boundary found in GOTHiC data between chromosome " + chrs[pos-1] + " and " + chrs[pos] + "\n")
		
	#	print ("position of chr: " + str(pos))
		#print ("position of chr: " + str(len(chrs)))

		rb = bSearchChr(chrs[pos], chrs[pos+1], 0, len(goth_data[0][:])-1)  
		if rb == -1:
			print("4")
                        exit("ELSE STATEMENT: No boundary found in GOTHiC data between chromosome " + chrs[pos] + " and " + chrs[pos+1] + "\n")

	#lb = lb + 2
	#rb = rb + 2

	#print ("lb: " + str(lb))
	#print ("rb: " + str(rb))

	#Now that we know the index boundaries of the chromosome, we need to find the bin that may or may not be relevant to chr1 interaction 



	#Here: for now, doing a basic linear search. May change this later if speed becomes an issue.	

	
		
	
	
		#Find second boundary

def buildGraph():
	#This function builds the interaction graph
	#First, create a vertex for each DM amplicon in the dmfinder_data dataframe
	dmID = 1
	#Here: grab a subset of the main dataframe
	#BELOW: THIS BELOW LINE DOES NOT WORK. MUST FIGURE OUT HOW TO GET A SUBSET OF DMFINDER INPUT DEPENDING ON THE DM INDEX VALUE 
	#dm = dmfinder_data.loc[dmfinder_data['index'] == "1"]

	global window_size #I don't like this. May have to make it read directly from data later on

	
#	robjects.r('dm <- dmfinder_data[dmfinder_data$index %in% c("dm' + str(dmID) + '")]')
	#Do a pairwise check on each DM amplicon to see if there is heightened activity

	global l
	global goth_data	
	global qthresh	
	global g	

	l = set(list(dmfinder_data[:][3]))

		
#	print ("dmfinder_data contents:")
	#print(dmfinder_data[4][0])
	#print(dmfinder_data[4][1])
	#print(dmfinder_data[4][2])
	#print(dmfinder_data[4][3])
	#print(dmfinder_data[4][4])
	for idx in l:
		#print ("idx: " + str(idx))
	
		for i in range(0, len(dmfinder_data[0][:])-1):
			#print ("In 1st for")
			for j in range(i + 1, len(dmfinder_data[0][:])):
				#print ("BEGIN DMFINDER DATA:")
				#print (dmfinder_data[i][0])
				#print (dmfinder_data[j][0])
				#print (dmfinder_data[i][0] + " " + dmfinder_data[j][0])
				#Here: we are checking if there is heightened interaction between all pairs of amplicons in a DM		
				if dmfinder_data[3][i] == idx and dmfinder_data[3][j] == idx:				
					#print ("In IF")
					#i is amplicon 1, j is amplicon 2
					#HERE: find the bin in the gothic data that amp1 belongs to
					#To find the bin, we first must find the location of the first and last record of the chromosome
					chromBinSearch(dmfinder_data[0][i])
					
					
					#print ("Length of goth_data: " + str(len(goth_data[1][:])))
					#print ("Length of dmfinder_data: " + str(len(dmfinder_data[1][:])))
					#print ("k: "  + str(k))
					#print ("i: " + str(i))
					#print ("j: " + str(j))
					#HERE: Boundaries in GOTHIC data found for given chromosome. Now do linear search to find the bin
					for k in range(lb, rb):
						#print ("in innermost FOR")
						if goth_data[1][k] < dmfinder_data[2][i] and goth_data[1][k] + window_size > dmfinder_data[1][i] and goth_data[3][k] < dmfinder_data[2][j] and goth_data[3][k] + window_size > dmfinder_data[1][j] and goth_data[10][k] <= qthresh and dmfinder_data[0][i] <= dmfinder_data[0][j] and goth_data[0][k] == dmfinder_data[0][i] and goth_data[2][k] == dmfinder_data[0][j]:
							#print ("Row is " + str(goth_data[0][k] + " " + str(goth_data[1][k]) + " " + str(goth_data[2][k]) + " " + str(goth_data[3][k]) + " " + str(goth_data[4][k]) + " " + str(goth_data[10][k])))
							#print ("qthresh 2: " + str(qthresh))
							#print("q-value: " + str(goth_data[10][k]))
							u = dmfinder_data[0][i] + "\n" + str(dmfinder_data[1][i]) + "\n" + str(dmfinder_data[2][i])
							v = dmfinder_data[0][j] + "\n" + str(dmfinder_data[1][j]) + "\n" + str(dmfinder_data[2][j])
							g.add_edge(u,v)
							#INTERACTION FOUND
							print ("HIC Interaction found between " + dmfinder_data[0][i] + ":" + str(dmfinder_data[1][i]) + "-" + str(dmfinder_data[2][i]) + " and " + dmfinder_data[0][j] + ":" + str(dmfinder_data[1][j]) + "-" + str(dmfinder_data[2][j]))
							break
					if dmfinder_data[0][i] > dmfinder_data[0][j]: #Same thing but now check for other interacting region if chr_i is lex. > chr_j. This is a quirk of how Gothic data is organized
							#Find bounding regions for other chromosome
						(lb_j, rb_j) = chromBinSearch2(dmfinder_data[0][j])	
						if(bSearchFirstBin(dmfinder_data[1][j], lb_j, rb_j) != -1):	#Search was successful
							#Put edgi
							for k in range(lb_j, rb_j):
		                                                #print ("in innermost FOR")
                		                                if goth_data[1][k] < dmfinder_data[2][j] and goth_data[1][k] + window_size > dmfinder_data[1][j] and goth_data[3][k] < dmfinder_data[2][i] and goth_data[3][k] + window_size > dmfinder_data[1][i] and goth_data[10][k] <= qthresh and goth_data[0][k] == dmfinder_data[0][i] and goth_data[2][k] == dmfinder_data[0][j]:
									#print ("qthresh 3: " + str(qthresh))
									#print("q-value: " + str(goth_data[10][k]))
									u = dmfinder_data[0][i] + "\n" + str(dmfinder_data[1][i]) + "\n" + str(dmfinder_data[2][i])
	                        	                                v = dmfinder_data[0][j] + "\n" + str(dmfinder_data[1][j]) + "\n" + str(dmfinder_data[2][j])
        		                                                g.add_edge(u,v)
									print ("HIC Interaction found between " + dmfinder_data[0][i] + ":" + str(dmfinder_data[1][i]) + "-" + str(dmfinder_data[2][i]) + " and " + dmfinder_data[0][j] + ":" + str(dmfinder_data[1][j]) + "-" + str(dmfinder_data[2][j]))
									
									break

					#Check if same dm
					#now check if heigtened interaction in gothic data
					
					#First do binary search to get correct chromosome
					
				else:
					break	#Don't test pairwise interaction if amplicon belongs to different predicted DM 
		
		plt.figure(figsize=(12,12))
		#plt.axis('off')
		pos = nx.spring_layout(g)
		nx.draw_networkx(g, pos, with_labels = True, font_size = 12)
		#nx.draw(g, with_labels = True)
		plt.savefig(myfile + "." + idx + ".png") # save as png
		#plt.show()
		#plt.show() # display	
	
		print ("Max clique graph file written to " + myfile + "." + idx + ".png")
		print ("Finding maximal cliques...")
		findCliques(g)
		g.clear()
		

def findCliques(g):
	global dmCounter
	line = ""

	if path.exists(myfile + ".disambig.bed"):
		os.remove(myfile + ".disambig.bed")

	f = open(myfile + ".disambig.bed", 'a+')
	#print ("Maximal cliques are: ")
	
	for D in list(find_cliques(g)):	
		for d in D:
			x = re.split('\n', d)
			#print ("d is : " + d)	
			#print ("x is : " + str(x))
			f.write(x[0] + "\t" + x[1] + "\t" + x[2] + "\t" + "DM" + str(dmCounter) + "\n")	
		dmCounter += 1
	#print (list(find_cliques(g)))
	f.close()

print ("Importing BED file...")
countDM_DMFinder() #DMFinder output is exported to this list
l = list (dmfinder_data[3][:])
indices = (set(l))

print ("Importing GOTHiC CSV file...\n")
importGothic()
buildGraph()

print("Output BED file is " + myfile + ".disambig.bed")
#print mylist[1]

#G.add_node([1,2])i
