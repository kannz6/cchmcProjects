#!/usr/bin/env python
import os
import re
import subprocess
import sys

from os import path

bamFiles = set()
currentDirectoryFiles = []; doneFileNames = [];

def joinPaths ( rt, fname ):
    filepath = os.path.join( rt, fname )
    doneFileNames.append( fname )  # Add it to the list.

def getFileNames( directory ):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    # Walk the tree.
    [ ([joinPaths(root, f) for f in files],[ currentDirectoryFiles.append(d) for d in directories]) for root, directories, files in os.walk( directory ) ]

    return doneFileNames  # Self-explanatory.

def lsfJobWrapper(w,M,n,R):
    return "#!/bin/bash\n#BSUB -W %s\n#BSUB -M %s\n#BSUB -n %s\n#BSUB -R \"span[%s]\"\n\n"%(w,M,n,R)

script = lsfJobWrapper("01:00","2000","1",("ptile=%s"%("1")))
directoryFileNames = getFileNames( os.getcwd() )
listOfBamFileNames = filter( (lambda x : re.match( r'.*aligned-sorted.*bam$', x) ), directoryFileNames )

with open("handleMultiKeyBams.sh", "w+") as handleMultiKeyBams:
    handleMultiKeyBams.write(script)
    handleMultiKeyBams.write("#list of bam files: {0}\n\n".format(listOfBamFileNames))

    for x in listOfBamFileNames:
        subDir = [ "{0}/{1}/{2}".format(os.getcwd(),i,x) for i in currentDirectoryFiles if os.path.isfile("{0}/{1}/{2}".format(os.getcwd(),i,x)) ]
        bamFiles.add(max(subDir,key=os.path.getsize))

    [ handleMultiKeyBams.write("mv {0}* {1}\n".format(bamFile,os.getcwd())) for bamFile in bamFiles ]
    handleMultiKeyBams.write("\necho \"Finished moving bam files!\" > {0}\n\nexit\n".format("completed-handleMultiKeyBams.txt"))
