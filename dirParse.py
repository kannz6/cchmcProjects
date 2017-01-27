# http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
import os
import re

file_paths = []  # List which will store all of the full filepaths.

def joinPaths ( rt, fname ):
	filepath = os.path.join( rt, fname )
	file_paths.append( filepath )  # Add it to the list.

def get_filepaths( directory ):
	"""
	This function will generate the file names in a directory 
	tree by walking the tree either top-down or bottom-up. For each 
	directory in the tree rooted at directory top (including top itself), 
	it yields a 3-tuple (dirpath, dirnames, filenames).
	"""
	# file_paths = [] # List which will store all of the full filepaths.

	# Walk the tree.
	[ [joinPaths(root, x) for x in files] for root, directories, files in os.walk( directory ) ]
	# for root, directories, files in os.walk( directory ):
	# 	[joinPaths(root, x) for x in files]

		# for filename in files:
		# 	# Join the two strings in order to form the full filepath.
		# 	filepath = os.path.join( root, filename )
		# 	file_paths.append( filepath )  # Add it to the list.

	return file_paths  # Self-explanatory.

# these two containers are duplicatedon lines 37 and 39 for scope. At this level, scope is all files for uniqueness
# at the next level, scope is local to the file being parsed
#container to store all the relation strings and FID's if needed
# kin_file_full_relation_list = [];  kin_file_fid_value_list = []; 
#  using these variables from refactor
kin_file_full_relation_list2 = []; kin_file_fid_value_list2 = [];

def appendFIDValueList( v ):
	if ( kin_file_fid_value_list2.count( v ) == 0 ):
		kin_file_fid_value_list2.append( v )
# ro is the relation string in its original order, rr is with FID1 and FID2 reveresed
def appendRelationStringList( ro, rr ):
	if ( kin_file_full_relation_list2.count( ro ) == 0 and kin_file_full_relation_list2.count( rr ) == 0 ):
		kin_file_full_relation_list2.append( ro )

def createRelationStrings( variablesList ):
	relation_str = variablesList[0] + " " + variablesList[2] + " " + variablesList[7]
	reverse_relation_str = variablesList[2] + " " + variablesList[0] + " " + variablesList[7]
	appendFIDValueList( variablesList[0] )
	appendFIDValueList( variablesList[2] )
	appendRelationStringList( relation_str, reverse_relation_str )

def parseKinFile( kinFile ):
	with open( kinFile ) as kf:
		for line in kf:
		#find relevant data in file
			kinFileVariables = line.split()
			if (kinFileVariables[0] == 'FID1'):
				continue

			createRelationStrings( kinFileVariables )

# Run the above function and store its results in a variable.   
full_file_paths = get_filepaths( "../../sandbox" )

# get all the files that match the file extension we are looking for
list_of_kin_files = filter( (lambda x : re.match( r'(.*\/.*kin0)', x) ), full_file_paths )
# using the list from the kin files, extract the variables we need
map( parseKinFile, list_of_kin_files )
#use this in the table
# for found_file in full_file_paths:
	#search for *.kin0 files
	# kinFileDir_regEx = re.match( r'(.*\/.*kin0)', found_file )
	# #container to store all the relation strings if needed
	# # kin_file_full_relation_list = [];
	# # #list of id's for sorting for unique values before entering into the table
	# # kin_file_fid_value_list = [];

	# if kinFileDir_regEx:

		# print kinFileDir_regEx.groups()

		# map( parseKinFile, kinFileDir_regEx.groups() )

		# for kin_file_in_directory in kinFileDir_regEx.groups():

		# 	with open( kin_file_in_directory ) as kin_file:

		# 		for line in kin_file:
		# 			#find relevant data in file
		# 			kinFileVariables = line.split()
		# 			if (kinFileVariables[0] == 'FID1'):
		# 				continue

		# 			createRelationStrings( kinFileVariables )					

	# in this scope, containers are local to the file only
	# for relationString in kin_file_full_relation_list:
	# 	print relationString

	# for fID in kin_file_fid_value_list:
	# 	print fID
# moved out for scope so that all files can be checked for uniqueness
# for relationString in kin_file_full_relation_list:
# 	print relationString

# for fID in kin_file_fid_value_list:
# 	print fID

print "printing kin file variables:\n"
for relationString in kin_file_full_relation_list2:
	print relationString

for fID in kin_file_fid_value_list2:
	print fID

