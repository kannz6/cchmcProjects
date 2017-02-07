from dirParse import dirParse
import re

# testing the class
tester = dirParse("../../sandbox")
tester.parseDirectory()
# tester.printRelationList()
# tester.printFIDValues()
# >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
# corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively.
# tester.printRelationDictionary()
# print "\n"
# "projects/pcgc/prod/data/expedat/public/Other/exome/fastq/"
tester.parseExomeFilePaths( "../../b2bscripts/TrioCheck3/parseFiles/test2.txt" )
tester.validateTrio()
tester.printBlindIdRelationData()
# for r in tester.getRelationDictionary():
# 	idVals = r[0].split() #r[0] contains the FID's/key in the dictionary
# 	relVals = r[1].split() #r[1] contains the relation/value in the dictionary
# 	print "_______________________________________"
# 	if ( float(relVals[2] ) > 0.354 ):
# 		print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: Duplicate/MZ Twin Relation" #+ r[1]
# 	elif ( float(relVals[2] ) <= 0.354 and float( relVals[2] ) >= 0.177 ):
# 		print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 1st-degree Relation" #+ r[1]
# 	elif ( float(relVals[2] ) < 0.177 and float( relVals[2] ) >= 0.0884 ):
# 		print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 2nd-degree Relation" #+ r[1]
# 	elif ( float(relVals[2] ) < 0.0884 and float( relVals[2] ) >= 0.0442 ):
# 		print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 3rd-degree Relation" #+ r[1]
# 	elif ( float( relVals[2] ) < 0.0442 ):
# 		print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 4th-degree Relation" #+ r[1]

# 	if ( float(relVals[2] ) > 0.0442 and ( ( idVals[1].endswith('-01') or idVals[1].endswith('-02') ) 
# 																	and ( idVals[4].endswith('-01') or idVals[4].endswith('-02') ) ) ):

		# print "Higher than expected value for parents relation!" 

