# import some required libraries
import sys
import os
import numpy as np
import re # for regexpression

# filenames (must be in same directory as python script) of the fastq files (not fastq.gz) to process
innames =  ["Undetermined_S0_L001_I1_001.fastq",
            "Undetermined_S0_L001_I2_001.fastq",
            "Undetermined_S0_L001_R1_001.fastq",
            "Undetermined_S0_L001_R2_001.fastq",
            "Undetermined_S0_L002_I1_001.fastq",
            "Undetermined_S0_L002_I2_001.fastq",
            "Undetermined_S0_L002_R1_001.fastq",
            "Undetermined_S0_L002_R2_001.fastq"]

firstFileForLane = ["Undetermined_S0_L001_I1_001.fastq","Undetermined_S0_L002_I1_001.fastq"] # these are the first files for the individual lanes - this is going to help us figure out what files we need to create our files that contain all mismatches and their matches.

# provide the index markers for the individual libraries here
ref_index = {"lib": ["lib64237", "lib64238", "lib64239", "lib64240"], 
             "index": ["AACCACGCAT", "CCCACCACAA", "AAGGGTTTAC", "AGCGCCTTGC"],
             "index2": ["TAACCTGAAT", "AAGCGGAGGT", "CGCGTGAGTA", "CACTTCGTAC"],}

combined_indeces = []

# combine provided indeces to index1+index2 as they are found in fastq files
for index1, index2 in zip(ref_index["index"], ref_index["index2"]):
    combined_indeces.append("+".join([index1,index2]))

print("You provided the following markers:") 
for ilib, imarker in zip(ref_index["lib"], combined_indeces):
    print(ilib + ": " + imarker)

# combine all the index1+index2 strings into one large string - this is used for matching
index_string = "--".join(combined_indeces)

# now define a function that matches a recovered pattern (from fastq file) to the combined index string 
def my_matcher(pattern, indeces, maxN = 4):  
    # maxN = maximum number of Ns per index string (2 per fragment)

    # First verify that there aren't too many N's in the string that is to be matched
    split_index = re.split("\+", pattern)

    if len(re.findall("N",split_index[0])) <= maxN and len(re.findall("N",split_index[1])) <= maxN: 
        # replace Ns with wildcard characters, + with \+ (+ is a regexp operator and won't get recognized correctly)
        pattern = pattern.replace("N",".")
        pattern = pattern.replace("+","\+")

        # match the pattern agains the string of all the indeces
        matches = re.findall(pattern,indeces)
        if len(matches) == 1:
            # if there is exactly one match, return that match
            match = matches[0]
        else:
            # if there are 0 or more than 1 possible matches, abort
            match = []

    else:
        # if there are too many Ns in the pattern, abort
        match = []

    return match

# this is pattern that we use to process the first line in a fastq fragment (total of 4 lines),
# the first line contains flowcellid, sample position and the index that was read out
# we are interested in that index - index can only be A,C,G,N,T and is exactly 10 nucleotides lont
markers_pattern = re.compile("(.*)([ACNGT]{10}\+[ACNGT]{10})")

# loop over all fastq files
for inname in innames:
    print("Processing " + inname)

    # if the file processed is the first file for a specific lane, we want to track what index fragments were detected and what they were changed to
    # the header of each fastq file entry is the same across the files for the lane
    # matched and unmatched indeces are written to csv files for each lane
    if inname in firstFileForLane:
        track = True
        matchFileName = "matches_lane" + str(firstFileForLane.index(inname)+1) + ".csv"
        nomatchFileName = "nomatches_lane" + str(firstFileForLane.index(inname)+1) + ".csv"
        matchFile = open(matchFileName, "w")
        nomatchFile = open(nomatchFileName, "w")
        print("First file for lane " + str(firstFileForLane.index(inname)+1) + ": writing matches and mismatches to .csv file" )
    else:
        track = False

    # names and paths for fastq files to write to
    outfiles = ["./" + lib + "/" + inname for lib in ref_index["lib"]]
    outfiles.append("./undetermined/" + inname)

    f = [open(ofnames, "w") for ofnames in outfiles]

    # loop over lines in file
    # as files are very large (easily 100 GB), it is not feasible to load the file
    # instead we are stepping through them line by line
    with open(inname, 'r') as fh:
        while True: # a workaround to keep looping until EOF
            line = fh.readline()
            if not line:
                # if that was just read in was empty (as in no lines available anymore), end the loop
                break

            # match the line agains the above defined pattern to grab the index1+index2 string from the fastq fragment
            x = markers_pattern.match(line.rstrip())

            # try to match fastq index1+index2 to one of the combinations you provided at the beginning of the file
            y = my_matcher(x[2],index_string)

            if y:
                # if there was a match, replace the fastq line and write it to the output file
                libidx = combined_indeces.index(y)
                NewLine = x[1] + y + "\n"
                f[libidx].write(NewLine)

                if track:
                    # if first file of lane, write fastq index1+index2 plus the index1+index2 it was replaced by to the csv file
                    matchFile.write(x[2] + ',' + y + "\n")

            else:
                # if no match, write to the undetermined file
                libidx = len(f)-1
                f[libidx].write(line)
                
                if track:
                    # if first file of lane, write fastq index1+index2 to the csv file that collects all the unmatchable strings
                    nomatchFile.write(x[2] + "\n")

            # copy over the next three lines (4 lines per fastq fragment per file)
            f[libidx].write(fh.readline())
            f[libidx].write(fh.readline())
            f[libidx].write(fh.readline())

    # properly close files after processing
    if track:
        matchFile.close()
        nomatchFile.close()
    for fi in f:
        fi.close()

print("All finished, congratz!")