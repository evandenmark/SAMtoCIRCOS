#!/usr/bin/python
###############################################################################
#
#    genomeCoverage.py version 1.0
#    
#    Calculates windowed coverage across a genome and finds mate
#    matches and positions for each genome scaffold using a SAM file
#
#    Calculates GC content across genome and gaps of unknown bases
#    for each genome scaffold using an optional FASTA file
#
#    Copyright (C) 2014 Matthew Neave & Evan Denmark
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

"""
NOTES ABOUT THIS PROGRAM

The following files are produced when only a SAM file is provided:

The coverage file produced gives the scaffold name, each subunit within each scaffold, and how many
times a read maps to that subunit. If no reads are mapped to the subunit, it will be given a value
of zero. In most cases (but not all) when there are consecutive subunits within a scaffold of 
coverage zero, there is a gap spanning these subunits.

The ends file produced gives the original scaffold with the original mates start and end position on
the original scaffold as well as its mate's scaffold start and position. For this version of the 
program, mates are 100 base pairs long because the program used to generate the mates (HiSeq) creates
mates of this length.

The karyotype file produced gives the scaffold name, its length, and the color in which the scaffold
will appear when it is run through a visualization software.



The following files are produced when the required SAM file and optional FASTA file are provided:

The GC file produced gives the scaffold name, each subunit within each scaffold, and the GC content of 
each scaffold. In most cases when the GC content of consecutive subunits is zero, it is due to a gap
spanning these scaffolds. In addition, if a windowSize is not specified in the command line, the default
is 1000. Therefore, one would expect every subunit to have a GC content percentage to no more than one 
decimal place (ex. 541 GC out of 1000 base pairs results in a 54.1% GC content). However, in many cases,
the GC content goes far beyond 1 decimal place because of gaps within the subunit. Some gaps may be only
a single nucleotide. Because this program does not count gaps as nucleotides when calculating GC content,
this results in a fraction with a denominator other than 1000, giving a percentage with many decimals.
Of course, this only applies when the default 1000 is used as the windowSize.

The gap file produced gives the scaffold in which the gap resides and the start and end position of the 
gap. This program defines a gap as unknown nucleotides, given as "N" in the FASTA file.



In addition to the file produced the following will be displayed to the user:

 - The number of reads in the SAM file that do not match to any scaffolds
 
 - Warnings if any scaffolds did not have any reads matched to them

 - Warnings if the scaffolds provided in your SAM are different than those in the FASTA (possibly due to a recombination of your scaffolds with outside software) 

"""

import argparse

parser = argparse.ArgumentParser(description = 'Please provide a SAM file (required), a windowSize (optional, default 1000), an end_size (optional, default 500), and a FASTA file (optional). ')
parser.add_argument('samfile',help= 'a sam file with your data')
parser.add_argument('-windowSize',default = 1000, type=int, help= 'window size for coverage subunits')
parser.add_argument('-end_size', default = 500, type=int, help = 'distance from end to calculate links')
parser.add_argument('-fasta', default = None, help= 'a fasta file to calculate gc content and gaps within the scaffolds')


samfile = parser.parse_args().samfile
windowSize = parser.parse_args().windowSize
end_size = parser.parse_args().end_size
fasta= parser.parse_args().fasta

# use argparse module to get input SAM file and size of windows

def genomeCoverage(samfile, windowSize, end_size, fasta):
    """
    The first object this function returns is a dictionary within a dictionary. The larger dictionary has the scaffold as a
    key and a dictionary as the value. The smaller dictionary has the subunit (within each scaffold) as the key and coverage
    value (an integer) as the value.
    The dictionary (scaffold_dict) is in the form {scaffold: {subunit:coverage}}.

    The second object this function returns is a dictionary with the scaffold as a key and
    a list containing (in this order) the map position of the original read and a list of the read's mate
    scaffold and the mate's position on the original scaffold.
    The dictionary (end_piece_dict) is in the form {scaffold:[original read position, [mate scaffold, mate position]]}

    The third object this function returns is a dictionary (len_dict) that simply shows the lengths for each scaffold.
    The dictionary (len_dict) is in the form {scaffold:length} 
    
    The fourth object this function returns is a dictionary (gc_dict) that shows the GC content of each subunit of 
    size windowSize of each scaffold.
    The dictionary is in the form {scaffold:{subunit:gc_content}}

    The fifth object this function returns is a dictionary (gap_dict) that shows the gaps of unknown base in each scaffold.
    The dictionary is in the form {scaffold: [[gap_start, gap_end], [gap_start, gap_end]]}

    """
    
    # open SAM file
    # initialize our dictionaries that we will use

    the_file = open(samfile)

    # the coverage dictionary
    scaffold_dict = {}      
    # scaffold length dictionary
    len_dict = {}
    # runs near the edge dictionary
    end_piece_dict = {}
    didNotMap = 0
    largest_coverage = {} 
    for each_line in the_file:
        """
        The first few lines of the sam file do not include information about the reads, but rather
        information about the scaffold
        """
        
        each_line = each_line.split()
        if each_line[0] == '@SQ':
            
            scaffold_name = each_line[1].split(':')
            
            name = str(scaffold_name[1]).lstrip('[').rstrip(']')
            name = name.strip()
            
            line_len = each_line[2].split(':')
            length = str(line_len[1]).lstrip('[').rstrip(']')
            
            length = length.strip()
            len_dict[name] = length

            
        if each_line[0][0] == '@':
            
            continue
        
        else:

            scaffold = each_line[2]
            mapPosition = int(each_line[3])
            possible_equal = each_line[6]
            mate_map = int(each_line[7])
	    include = True
               
            """
            if the original read is at the end of a scaffold and it has a mate on another scaffold, we want to remember
	    its position and mate, so we create the end_piece dictionary which, for each scaffold as a key, creates a 
	    giant list of smaller list for each read. Within these smaller lists, the first object is the map position 
	    of the original read and the second object is a list of the name of the mate scaffold and position of mate

            End_piece_dict dictionary is in the form {scaffold: [[mapPosition,[mate scaffold, mate position]], [mapPosition,[mate scaffold, mate position]] etc...]
        
            """
            if scaffold == '*':
                didNotMap += 1
                continue
                
            else:
		if possible_equal != '=':
			# the following ensures that the reciprocal end link is not included 
			if possible_equal in end_piece_dict:
				for each_mate in end_piece_dict[possible_equal]:
					if each_mate[0] == mate_map:
						include = False
						
		if include == True:
						
		        if end_piece(mapPosition,int(len_dict[scaffold]), possible_equal, end_size) and is_mate_at_end(int(len_dict[possible_equal]),mate_map, end_size):
	       			if scaffold not in end_piece_dict:
		               		end_piece_dict[scaffold] = []
			    	read_list = [mapPosition, [possible_equal, mate_map]]
		               	end_piece_dict[scaffold].append(read_list)
				
 	    
	    # coveragePosition is the name we will give each window of size windowSize (ex. if there are 10,000 bp and we have windowSize 1000, there will be 10 coverage positions of names 0-9)
            coveragePosition = mapPosition / windowSize
            coveragePosition = int(coveragePosition)
	    	    
	    if scaffold not in largest_coverage or coveragePosition > largest_coverage[scaffold]:
		largest_coverage[scaffold] = coveragePosition 

            if scaffold not in scaffold_dict:
                scaffold_dict[scaffold] = {}
      	    
            if coveragePosition in scaffold_dict[scaffold]:
                scaffold_dict[scaffold][coveragePosition] += 1
	        
            else:
                scaffold_dict[scaffold][coveragePosition]  = 1
		
    	    

    for each_scaffold in largest_coverage:
	for i in xrange(largest_coverage[each_scaffold]):
		if i not in scaffold_dict[each_scaffold]:
			scaffold_dict[each_scaffold][i] =0
			    

    print
    print 'For reference,', didNotMap, 'reads did not map to any scaffold.' 
    print 

    for each_scaffold in len_dict:
	if each_scaffold not in scaffold_dict:
		print
		print "WARNING! No reads mapped to", each_scaffold, "and therefore it will not be in your coverage file."
		print
    # A fasta file must be provided in order to make gc_content dictionary and gap dictionary
    if fasta != None:
    	fasta= open(fasta)
    	gc_dict = {}
	gap_dict = {}
    	for each_line in fasta:
		
    		if each_line[0] == '>':
			#name line
			each_line = each_line.lstrip('>')
			each_line = each_line.rstrip('\n')
			# HERE 'each_line' MUST BE THE NAME OF A SCAFFOLD
			# Some fasta files may differ in the way the original program gives the '>' line a name, so it may require alteration
			gc_dict[each_line] = {}
			name_line = each_line
			if each_line in scaffold_dict:				
				if name_line not in gap_dict:
					gap_dict[name_line] = []
				
			count = 0
			num_gc = 0.0
			num_actg = 0.0
			current_gap = False
			gc_content = 0.0	
		else:
			#sequence line
			each_line = each_line.rstrip('\n')
			for each_base in each_line:
				count +=1
				each_base = str(each_base)
				
				if current_gap == True:
					if each_base.upper() == 'N':
						gap_end+=1
						
					else:
						the_gap = [gap_start, gap_end]
						gap_dict[name_line].append(the_gap)
						current_gap = False
				elif each_base.upper() == 'N':
					gap_start = count
					gap_end = gap_start
					current_gap = True				

				else:
					if each_base.upper() == 'C' or each_base.upper() == 'G':
						num_gc +=1
						num_actg +=1
						gc_content = (float(num_gc)/float(num_actg))*100.0
					
					if each_base.upper() == 'A' or each_base.upper() == 'T':
						num_actg+=1
						gc_content = (float(num_gc)/float(num_actg))*100.0	
				if count%windowSize == 0:
					current_window=((int(count/windowSize))-1)
					gc_dict[name_line][current_window] = gc_content
					gc_content = 0.0
					num_gc = 0.0
					num_actg = 0.0
				elif count == int(len_dict[name_line]):
					gc_dict[name_line][current_window+1] = gc_content
					gc_content = 0.0
					num_gc = 0.0
					num_actg = 0.0	 
	
	for each_scaffold in scaffold_dict:
		if int((int(len_dict[each_scaffold]))/windowSize)+1 != len(gc_dict[each_scaffold]):
			print 
			print "WARNING! The scaffolds in the SAM file have different lengths than the scaffolds in the FASTA file."
			print
			break	  
	return scaffold_dict, end_piece_dict, len_dict, gc_dict, gap_dict
    else:
    	return scaffold_dict, end_piece_dict, len_dict


def end_piece(start_pos, length, possible_equal, end_size):
    """
    Determines if your read is at the end of a scaffold.

    """
    
    if (((length - start_pos) < end_size) or (start_pos < end_size)) and possible_equal != '=':
        #it is at the end and mate is on another scaffold
        return True
    else:
        #either the read is not near the end or the mate is on the same scaffold or both
        return False


def is_mate_at_end(length_mate_scaffold, position_of_mate, end_size):
    """
    Determines if the mate of your original read is near the end of its scaffold
    """

    if ((length_mate_scaffold - position_of_mate) < end_size) or (position_of_mate < end_size):
        #the mate is near the end of its scaffold
        return True
    else:
        return False


def make_files(samfile, windowSize, end_size, fasta):
    """
    Takes the organized data generated by genomeCoverage() and creates 3 files: a coverage file,
    a karyotype file, and an end file.
    """
    

    
    genome = genomeCoverage(samfile, windowSize, end_size, fasta)
    if fasta != None:
    	gap_dict = genome[4]
    	gc_dict = genome[3]
    len_dict = genome[2]
    end_piece_dict = genome[1]
    scaffold_dict = genome[0]
    
    scaff_list = []
    len_list = []
    end_list = []
    gc_list = []
    gap_list = []

    for scaffold in scaffold_dict:
        scaff_list.append(scaffold)
    scaff_list = sorted(scaff_list)

    for each_scaffold in len_dict:
        len_list.append(each_scaffold)
    len_list = sorted(len_list)

    for the_scaffold in end_piece_dict:
        end_list.append(the_scaffold)
    end_list = sorted(end_list)

    if fasta != None:
        for scaffolds in gc_dict:
            gc_list.append(scaffolds)
        gc_list= sorted(gc_list)

        for a_scaffold in in gap_dict:
            gap_list.append(a_scaffold)
        gap_list = sorted(gap_list)


    # MAKE THE COVERAGE FILE    
    new_file = open('output.coverage.txt', 'w')

    end=0
    for scaffold in scaff_list:
        length = int(len_dict[scaffold])
        for subunit in scaffold_dict[scaffold]:
            start = subunit*windowSize
            end = start + windowSize
            if end >= length:
                end = length
            coverage = scaffold_dict[scaffold][subunit]
            string = str(str(scaffold)+'\t' + str(start)+ '\t' +str(end)+ '\t' +str(coverage)+'\n')
            new_file.write(string)
    new_file.close()

    # MAKE THE KARYOTYPE FILE
    new_file2 = open('output.karyotype.txt', 'w')

    n=0
    for scaffold in len_list:
	if n == 7:
		n=0
        length = int(len_dict[scaffold])
        color = ['set1-7-qual-1','set1-7-qual-2','set1-7-qual-3','set1-7-qual-4','set1-7-qual-5','set1-7-qual-6','set1-7-qual-7']
        line = ('chr -' + '\t' + str(scaffold)+ '\t' +str(scaffold) +'\t' + str(0) + '\t' + str(length) + '\t' +str(color[n]) + '\n')
        new_file2.write(line)
	n+=1
    new_file2.close()

    # MAKE THE ENDS FILE

    new_file3 = open('output.ends.txt', 'w')
    for scaffold in end_list:
	for each_end in end_piece_dict[scaffold]:        
		original = str(scaffold)
        	original_start = str(each_end[0])
        	original_end = int(original_start) + 100
        	if original_end > int(len_dict[scaffold]):
            		original_end = int(len_dict[scaffold])
		original_end = str(original_end)
        	mate = str(each_end[1][0])
        	mate_start = str(each_end[1][1])
        	mate_end = str(int(mate_start)+100)
        	if int(mate_end) > int(len_dict[mate]):
            		mate_end = str(len_dict[mate])

        	output = str(original + '\t' + original_start + '\t' + original_end + '\t' + mate + '\t' + mate_start + '\t' +mate_end + '\n')

        	new_file3.write(output)
    new_file3.close()


    if fasta != None:

	    # MAKE THE GC FILE
	    new_file4 = open('output.gc.txt', 'w')
	    for scaffold in gc_list:
		the_scaffold = str(scaffold)
		length = int(len_dict[scaffold])
	    	for subunit in gc_dict[scaffold]:
			start = int(subunit)*windowSize
			end = start+windowSize
			if end > length:
				end = length
			end = str(end)
			start= str(start)
			gc_content = str(gc_dict[scaffold][subunit])
			output= str(the_scaffold + '\t' + start+ '\t'+ end + '\t' + gc_content + '\n')
			new_file4.write(output)
	    new_file4.close()
			

	    # MAKE THE GAP FILE
	    new_file5 = open('output.gaps.txt', 'w')
	    for scaffold in gap_list:
		the_scaffold = str(scaffold)
		for gap in gap_dict[scaffold]:
			gap_start=str(gap[0])
			gap_end=str(gap[1])
			output= str(the_scaffold + '\t' + gap_start + '\t' + gap_end + '\n')
			new_file5.write(output)
	    new_file5.close()
    
make_files(samfile, windowSize, end_size, fasta)
