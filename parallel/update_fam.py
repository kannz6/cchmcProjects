#!/usr/bin/env python
#Andrew Rupert 4-6-17 
import fileinput
import os
import re
import sys
######################
#bulk vcf (yale)
######################
# r_parse_base_id=re.compile(r'^(.*?)(-\d\d)?$')
# new_columns = []
path = "tmp_fam.txt"
# tmpFamFileWriter = open("tmp_fam.txt", "w+")
# for line in fileinput.input():
#     line=line.strip()
#     columns=line.split(' ')
#     new_columns=columns[:]
#     #print(columns)
#     #print(new_columns)
#     match = r_parse_base_id.match(columns[0])

#     family_id = match.group(1)
#     gender = '0'
#     mother_id = '0'
#     father_id = '0'
#     if match.group(2) == '-01':
#         gender = '2'
#     elif match.group(2) == '-02':
#         gender = '1'
#     elif match.group(2) is None:
#         mother_id = family_id + '-01'
#         father_id = family_id + '-02'
#     else:
#         raise ValueError('HEEELP')

#     new_columns[0] = family_id
#     new_columns[2] = mother_id
#     new_columns[3] = father_id
#     new_columns[4] = gender
#     tmpFamFileWriter.write(' '.join(new_columns))
#     tmpFamFileWriter.write("\n")
#     # print(' '.join(new_columns))#bulk yale vcf use > to write to file in
# tmpFamFileWriter.close()
######################
doneFileNames = [];
curretDirectoryFiles = [];
def joinPaths ( rt, fname ):
    filepath = os.path.join( rt, fname )
    doneFileNames .append( fname )  # Add it to the list.

def getFileNames( directory ):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    # Walk the tree.
    [ ([joinPaths(root, f) for f in files],[ curretDirectoryFiles.append(d) for d in directories]) for root, directories, files in os.walk( directory ) ]

    return doneFileNames  # Self-explanatory.
######################
#fastq to kin pipeline
######################
# r_parse_base_id=re.compile(r'^(.*?)(-\d\d)?$')
r_parse_base_id=re.compile(r'^((.*?)(-\d+|(-\d+-\d)))?$')
new_columns = []
path = "tmp_fam.txt"
tmpFamFileWriter = open("tmp_fam.txt", "w+")
for i,line in enumerate(fileinput.input()):
    line=line.strip()
    columns=line.split(' ')
    new_columns=columns[:]
    # print(columns)
    # print(new_columns)
    _id = new_columns[1].split(".")
   # print(_id)
    _id = _id[1]
    family_id = _id; filesDict = {};
    match = r_parse_base_id.match(columns[0])
    family_id = match.group(1)
    gender = '0'
    mother_id = '0'
    father_id = '0'
    # print("family id: {0}\nid: {1}\nmatch group 0: {2}\n".format(family_id,_id,match.group(0)))
    if (i == 0 and _id.count("-") == 2):

        directoryFileNames = getFileNames( os.getcwd() )
        listOfBamFileNames = filter( (lambda x : re.match( r'aligned-sorted.*bam$', x) ), directoryFileNames )

        [ filesDict.update({x.split(".")[1][0:7]:x.split(".")[1]+":"+x}) if x.split(".")[1].count("-") == 2  else filesDict.update({x.split(".")[1][0:10]:x.split(".")[1]+":"+x}) for x in listOfBamFileNames if ( len(x.split(".")[1]) == 9 and os.path.getsize(x) > filesDict.get(x.split(".")[1][0:7]) ) or ( len(x.split(".")[1]) == 12 and os.path.getsize(x) > filesDict.get(x.split(".")[1][0:10]))]
        # sys.exit("filesDict: {0}".format(filesDict))
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[2] = filesDict.get(family_id+"-01").split(":")[0]
        new_columns[3] = filesDict.get(family_id+"-02").split(":")[0]
        new_columns[4] = gender

    elif(i == 0 and _id.count("-") == 1):
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[2] = family_id + '-01'
        new_columns[3] = family_id + '-02'
        new_columns[4] = gender
    elif (i % 2) == 1:
        if ("01" in _id.split("-")):
            gender = '2'
        elif ("02" in _id.split("-")):
            gender = '1'
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[4] = gender
    elif (i % 2 ) == 0:
        if ("01" in _id.split("-")):
            gender = '2'
        elif ("02" in _id.split("-")):
            gender = '1'
        new_columns[0] = family_id
        new_columns[1] = _id
        new_columns[4] = gender
    tmpFamFileWriter.write(' '.join(new_columns))
    tmpFamFileWriter.write("\n")
    # print(' '.join(new_columns))#bulk yale vcf use > to write to file in
tmpFamFileWriter.close()
######################
