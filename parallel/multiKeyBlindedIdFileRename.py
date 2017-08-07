#!/usr/bin/env python
import os
import re
import subprocess

#################################################################################
#8-2-17
doneFileNames = [];
curretDirectoryFiles = [];
filesDict = {};
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
    [ ([joinPaths(root, f) for f in files],[ curretDirectoryFiles.append(d) for d in directories]) for root, directories, files in os.walk( directory ) ]

    return doneFileNames  # Self-explanatory.

def lsfJobWrapper(w,M,n,R):
    return "#!/bin/bash\n#BSUB -W %s\n#BSUB -M %s\n#BSUB -n %s\n#BSUB -R \"span[%s]\"\n\n"%(w,M,n,R)

directoryFileNames = getFileNames( os.getcwd() )
listOfBamFileNames = filter( (lambda x : re.match( r'.*aligned-sorted.*bam$', x) ), directoryFileNames )
with open("renameBam.sh", "w+") as renameBam:
    #stopped here 8-4-17 need to handle the
    [ filesDict.update({x.split(".")[1][0:7]:x.split(".")[1]+":"+x}) if x.split(".")[1].count("-") == 2  else filesDict.update({x.split(".")[1][0:10]:x.split(".")[1]+":"+x}) for x in listOfBamFileNames if ( len(x.split(".")[1]) == 9 and os.path.getsize(x) > filesDict.get(x.split(".")[1][0:7]) ) or ( len(x.split(".")[1]) == 12 and os.path.getsize(x) > filesDict.get(x.split(".")[1][0:10]))]
    # print("directoryFileNames: {0}\nlistOfBamFileNames: {1}\nfilesDict: {2}".format(directoryFileNames,listOfBamFileNames,filesDict))
    script = lsfJobWrapper("01:00","2000","1",("ptile=%s"%("1")))
    renameBam.write(script)
    # [ renameBam.write("mv {0} {1}\n".format(str(v.split(":")[1]) ,"aligned-sorted.{0}-test.bam".format(x))) for x,v in filesDict if ( len(x) == 9 and v.split(":")[0].count("-") == 2 ) or ( len(x) == 12 and v.split(":")[0].count("-") == 3 ) ]
    [ renameBam.write("mv {0} {1}\n".format(str(v.split(":")[1]) ,"aligned-sorted.{0}.bam".format(x))) for x,v in filesDict.iteritems() ]
    renameBam.write("\n\nexit\n")

