# http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
import os
import re

def get_filepaths( directory ):
    """
    This function will generate the file names in a directory 
    tree by walking the tree either top-down or bottom-up. For each 
    directory in the tree rooted at directory top (including top itself), 
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk( directory ):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join( root, filename )
            file_paths.append( filepath )  # Add it to the list.

    return file_paths  # Self-explanatory.

# these two containers are duplicatedon lines 37 and 39 for scope. At this level, scope is all files for uniqueness
# at the next level, scope is local to the file being parsed
#container to store all the relation strings and FID's if needed
kin_file_full_relation_list = []; kin_file_fid_value_list = [];

# Run the above function and store its results in a variable.   
full_file_paths = get_filepaths( "../../sandbox" )

#use this in the table
for found_file in full_file_paths:
	#search for *.kin0 files
	kinFileDir_regEx = re.match( r'(.*\/.*kin0)', found_file )
	#container to store all the relation strings if needed
	# kin_file_full_relation_list = [];
	# #list of id's for sorting for unique values before entering into the table
	# kin_file_fid_value_list = [];

	if kinFileDir_regEx:

		for kin_file_in_directory in kinFileDir_regEx.groups():

			#show which directory we're in
			with open( kin_file_in_directory ) as kin_file:
				for line in kin_file:
					#find relevant data in file
					kinFileVariables = line.split()
					if (kinFileVariables[0] == 'FID1'):
						continue

					relation_str = kinFileVariables[0] + " " + kinFileVariables[2] + " " + kinFileVariables[7]
					reverse_relation_str = kinFileVariables[2] + " " + kinFileVariables[0] + " " + kinFileVariables[7]
					#check if the FID's already exist in the list, if not...append them.
					if ( kin_file_fid_value_list.count( kinFileVariables[0] ) == 0 ):
						kin_file_fid_value_list.append( kinFileVariables[0] )
					if ( kin_file_fid_value_list.count( kinFileVariables[2] ) == 0 ):
						kin_file_fid_value_list.append( kinFileVariables[2] )

					#check if the relation exists in the list already, if not, append it. Make sure the reverse has not already been added
					if ( kin_file_full_relation_list.count( relation_str ) == 0 and kin_file_full_relation_list.count( reverse_relation_str ) == 0 ):
						kin_file_full_relation_list.append( relation_str )
					# print "Original Order: " + relation_str + "\nReverse Order: " + reverse_relation_str
	# in this scope, containers are local to the file only
	# for relationString in kin_file_full_relation_list:
	# 	print relationString

	# for fID in kin_file_fid_value_list:
	# 	print fID
# moved out for scope so that all files can be checked for uniqueness
for relationString in kin_file_full_relation_list:
	print relationString

for fID in kin_file_fid_value_list:
	print fID

