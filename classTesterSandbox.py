from dirParse import dirParse
import re

# testing the class
tester = dirParse("../../sandbox")
tester.parseDirectory()
# tester.printRelationList()
# tester.printFIDValues()
# >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] 
# corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively. 
for r in tester.getRelationList():
	relVals = r.split()
	if ( float(relVals[2] ) > 0.354 ):
		print "Duplicate/MZ Twin Relation: " + r
	if ( float(relVals[2] ) <= 0.354 and float( relVals[2] ) >= 0.177):
		print "1st-degree Relation: " + r
	if ( float(relVals[2] ) < 0.177 and float( relVals[2] ) >= 0.0884):
		print "2nd-degree Relation: " + r
	if ( float(relVals[2] ) < 0.0884 and float( relVals[2] ) >= 0.0442):
		print "3d-degree Relation: " + r