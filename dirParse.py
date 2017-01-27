# http://stackoverflow.com/questions/3207219/how-to-list-all-files-of-a-directory-in-python
import os
import re

class dirParse:

	def __init__(self, directory):
		self.dir = directory
		self.file_paths = []  
		self.kin_file_full_relation_list = []; 
		self.kin_file_fid_value_list = [];
		self.kin_file_relation_dict = {}

	def joinPaths ( self, rt, fname ):
		filepath = os.path.join( rt, fname )
		self.file_paths.append( filepath )  # Add it to the list.

	def get_filepaths( self, directory ):
		"""
		This function will generate the file names in a directory 
		tree by walking the tree either top-down or bottom-up. For each 
		directory in the tree rooted at directory top (including top itself), 
		it yields a 3-tuple (dirpath, dirnames, filenames).
		"""
		# Walk the tree.
		[ [self.joinPaths(root, f) for f in files] for root, directories, files in os.walk( directory ) ]

		return self.file_paths  # Self-explanatory.

	def appendFIDValueList(self, v ):
		if ( self.kin_file_fid_value_list.count( v ) == 0 ):
			self.kin_file_fid_value_list.append( v )
	# ro is the relation string in its original order, rr is with FID1 and FID2 reveresed
	def appendRelationStringList( self, ro, rr ):
		if ( self.kin_file_full_relation_list.count( ro ) == 0 and self.kin_file_full_relation_list.count( rr ) == 0 ):
			self.kin_file_full_relation_list.append( ro )

	def addRelationToDictionary(self, k1, k2, v):
		if( not self.kin_file_relation_dict.has_key( k1 ) ):
			self.kin_file_relation_dict.update( { k1 : v } )
		elif( not self.kin_file_relation_dict.has_key( k2 ) ):
			self.kin_file_relation_dict.update( { k2 : v } )

	def createRelationStrings( self, variablesList ):
		recordID_Part1  = ""
		recordID_Part2 = ""
		bothIDsFound = True
		# filter( (lambda x : re.match( r'.*-([0-9]+)-.*', x) ), full_file_paths )
		# aligned-sorted.1-13996-01.0.bam aligned-sorted.1-13996-02.1.bam 0.0907
		# aligned-sorted.1-13996-01.0.bam aligned-sorted.1-13996.2.bam 0.1732
		# aligned-sorted.1-13996-02.1.bam aligned-sorted.1-13996.2.bam 0.1723
		recordID_regExPart1 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*).*', variablesList[0] ) 
		recordID_regExPart2 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*|[0-9]+-[0-9]+)\..*', variablesList[2] ) 
		if recordID_regExPart1:
			recordID_Part1 += recordID_regExPart1.group(1)
			# print "id Part 1: %s"% recordID_Part1
		else:
			bothIDsFound = False

		if recordID_regExPart2:
			recordID_Part2 += recordID_regExPart2.group(1)
			# print "id Part 2: %s"% recordID_Part2
		else:
			bothIDsFound = False

		relation_str = variablesList[0] + " " + variablesList[2] + " " + variablesList[7]
		reverse_relation_str = variablesList[2] + " " + variablesList[0] + " " + variablesList[7]
		self.appendFIDValueList( variablesList[0] )
		self.appendFIDValueList( variablesList[2] )
		self.appendRelationStringList( relation_str, reverse_relation_str )
		# only add if the id's have been found so not to throw an error when calling method
		if ( bothIDsFound ):
			self.addRelationToDictionary( recordID_Part1, recordID_Part2, relation_str )

	def parseKinFile( self, kinFile ):
		with open( kinFile ) as kf:
			for line in kf:
			#find relevant data in file
				kinFileVariables = line.split()
				if (kinFileVariables[0] == 'FID1'):
					continue

				self.createRelationStrings( kinFileVariables )

	def printRelationList( self ):
		print "\nRelation List:\n___________________________"
		for r in self.kin_file_full_relation_list:
			print r
		print "___________________________\nEnd of Relation List"

	def printRelationDictionary( self ):
		print "\nRelation Dictionary:\n___________________________"
		# for k in self.kin_file_relation_dict.keys():
		# 	print "key: %s, value: %s" % (k, self.kin_file_relation_dict.get(k))
		dictList = self.kin_file_relation_dict.items()
		dictList.sort()
		for i in dictList:
			print "key: %s, value %s"% (i[0], i[1])

		print "___________________________\nEnd of Relation Dictionary"

	def getRelationList( self ):
		return self.kin_file_full_relation_list

	def printFIDValues( self ):
		print "\nFID Values:\n___________________________\n"
		for fid in self.kin_file_fid_value_list:
			print fid
		print "___________________________\nEnd of FID Values"

	def getFIDValues( self ):
		return self.kin_file_fid_value_list

	def parseDirectory( self ):
		# Run the above function and store its results in a variable.   
		full_file_paths = self.get_filepaths( self.dir )
		# get all the files that match the file extension we are looking for
		list_of_kin_files = filter( (lambda x : re.match( r'(.*\/.*kin0)', x) ), full_file_paths )
		# using the list from the kin files, extract the variables we need
		map( self.parseKinFile, list_of_kin_files )
