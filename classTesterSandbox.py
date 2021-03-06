from dirParse import dirParse
import re

# testing the class
tester = dirParse("../../sandbox")
tester = dirParse("../../sandbox/python/kinofiles/3-23-17")
# tester = dirParse("../../b2bscripts/TrioCheck3/mw_validation2_17_17")
# tester = dirParse("/scratch/kannz6/")
tester.parseDirectory()
# tester.printRelationList()
# tester.printFIDValues()
# >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
# corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively.
# tester.printRelationDictionary()
# tester.getRelationDictionary()
# print "\n"
# "projects/pcgc/prod/data/expedat/public/Other/exome/fastq/"
# tester.parseExomeFilePaths( "../../b2bscripts/TrioCheck3/parseFiles/test2.txt" )
# tester.parseExomeFilePaths()
tester.validateTrio()
print "\n"
# tester.printBlindIdRelationData()
print tester.getSpreadSheetFortmatedBlindIdRelationData()
