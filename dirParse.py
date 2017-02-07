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
		#data structures for creating the spreadsheet
		self.blindIdFilePathDict = {}
		self.allChildBlindIds = [];
		self.trioValidationDict = {}


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

	def appendChildFIDValueList(self, v):
		if( self.allChildBlindIds.count( v ) == 0 and not v.endswith("-01") and not v.endswith("-02") ):
			# print "Adding fid: %s to list of fid's"% v
			self.allChildBlindIds.append( v )

	def appendFIDValueList(self, v ):
		if ( self.kin_file_fid_value_list.count( v ) == 0 ):
			self.kin_file_fid_value_list.append( v )

	# ro is the relation string in its original order, rr is with FID1 and FID2 reveresed
	def appendRelationStringList( self, ro, rr ):
		if ( self.kin_file_full_relation_list.count( ro ) == 0 and self.kin_file_full_relation_list.count( rr ) == 0 ):
			self.kin_file_full_relation_list.append( ro )

	def addRelationToDictionary(self, k1, k2, v):
		dual_k1 = "FID1: " + k1 + " | FID2: " + k2; dual_k2 = "FID1: " + k2 + " | FID2: " + k1
		# if( not self.kin_file_relation_dict.has_key( k1 ) ):
		# 	self.kin_file_relation_dict.update( { k1 : v } )
		# elif( not self.kin_file_relation_dict.has_key( k2 ) ):
		# 	self.kin_file_relation_dict.update( { k2 : v } )
		if( not self.kin_file_relation_dict.has_key( dual_k1 ) ):
			self.kin_file_relation_dict.update( { dual_k1 : v } )
		elif( not self.kin_file_relation_dict.has_key( dual_k2 ) ):
			self.kin_file_relation_dict.update( { dual_k2 : v } )
		else:
			print "K->V already exists for: " + dual_k1 + ", and for: " + dual_k2 + "\nRelation: " + v

	def createRelationStrings( self, variablesList ):
		recordID_Part1 = "" ; recordID_Part2 = ""; recordID_Part3 = "" ; recordID_Part4 = ""
		bothIDsFound_1_2 = True; bothIDsFound_3_4 = True
		recordID_regExPart1 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*).*', variablesList[0] );
		recordID_regExPart2 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*|[0-9]+-[0-9]+)\..*', variablesList[2] ) 
		recordID_regExPart3 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[0] ); recordID_regExPart4 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[2] ) 

		if recordID_regExPart1:
			recordID_Part1 += recordID_regExPart1.group(1);	bothIDsFound_3_4 = False
						#append child id's only #2-7-17
			self.appendChildFIDValueList( recordID_regExPart1.group(1) )
			# print "id Part 1: %s"% recordID_Part1
		elif recordID_regExPart3:
			recordID_Part3 += recordID_regExPart3.group(1);	bothIDsFound_1_2 = False
			self.appendChildFIDValueList( recordID_regExPart3.group(1) )
		else:
			bothIDsFound_1_2 = False; bothIDsFound_3_4 = False


		if recordID_regExPart2:
			recordID_Part2 += recordID_regExPart2.group(1);	bothIDsFound_3_4 = False
			# print "id Part 2: %s"% recordID_Part2
		elif recordID_regExPart4:
			recordID_Part4 += recordID_regExPart4.group(1);	bothIDsFound_1_2 = False
		else:
			bothIDsFound_1_2 = False; bothIDsFound_3_4 = False

		relation_str = variablesList[0] + " " + variablesList[2] + " " + variablesList[7]
		reverse_relation_str = variablesList[2] + " " + variablesList[0] + " " + variablesList[7]
		self.appendFIDValueList( variablesList[0] ); self.appendFIDValueList( variablesList[2] )

		self.appendRelationStringList( relation_str, reverse_relation_str )
		# only add if the id's have been found so not to throw an error when calling method
		# if ( bothIDsFound ):
		# 	self.addRelationToDictionary( recordID_Part1, recordID_Part2, relation_str )
		if ( bothIDsFound_1_2 ):
			self.addRelationToDictionary( recordID_Part1, recordID_Part2, relation_str )
		elif ( bothIDsFound_3_4 ):
			self.addRelationToDictionary( recordID_Part3, recordID_Part4, relation_str )

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

	def getRelationDictionary( self ):
		dictList = self.kin_file_relation_dict.items()
		dictList.sort()
		return dictList

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

	def parseExomeFilePaths( self, file ):
		with open( file, "r" ) as exomeFileReaderO:
			exomeFileContent = exomeFileReaderO.readlines()
			exomeFileReaderO.close()
			exomeFileContent = [line.strip() for line in exomeFileContent] 
			
			for exomeFilePath in exomeFileContent:
				# print exomeFilePath
				blindIDPath_regEx = re.search( r'.*(1-[0-9]+-[0-9]+|1-[0-9]+).*', exomeFilePath );
				if blindIDPath_regEx:
					# print "Blind ID: %s"% blindIDPath_regEx.group(1)
					self.blindIdFilePathDict.update( { blindIDPath_regEx.group(1) : exomeFilePath } )
	
	def validateTrio( self ):
		for r in self.getRelationDictionary():
			idVals = r[0].split() #r[0] contains the FID's/key in the dictionary
			relVals = r[1].split() #r[1] contains the relation/value in the dictionary
			print "_______________________________________"
			if ( float(relVals[2] ) > 0.354 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: Duplicate/MZ Twin Relation" #+ r[1]
			elif ( float(relVals[2] ) <= 0.354 and float( relVals[2] ) >= 0.177 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 1st-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					self.trioValidationDict.update( { idVals[1] : "yes" } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					self.trioValidationDict.update( { idVals[4] : "yes" } )
			elif ( float(relVals[2] ) < 0.177 and float( relVals[2] ) >= 0.0884 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 2nd-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					self.trioValidationDict.update( { idVals[1] : "no" } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					self.trioValidationDict.update( { idVals[4] : "no" } )
			elif ( float(relVals[2] ) < 0.0884 and float( relVals[2] ) >= 0.0442 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 3rd-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					self.trioValidationDict.update( { idVals[1] : "no" } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					self.trioValidationDict.update( { idVals[4] : "no" } )
			elif ( float( relVals[2] ) < 0.0442 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 4th-degree Relation" #+ r[1]

			if ( float(relVals[2] ) > 0.0442 and ( ( idVals[1].endswith('-01') or idVals[1].endswith('-02') ) 
																			and ( idVals[4].endswith('-01') or idVals[4].endswith('-02') ) ) ):
				self.trioValidationDict.update( { idVals[1][:-3] : "no" } )

	def printBlindIdRelationData( self ):
		for blindID in self.allChildBlindIds:
			print "Child-id: %s, Child file: %s, Mom-id: %s, Mom file: %s, Dad-id: %s, Dad-file %s, Trio Validated: %s"%(blindID, self.blindIdFilePathDict.get(blindID), 
																							blindID + "-01", self.blindIdFilePathDict.get(blindID+"-01"),
																							blindID + "-02", self.blindIdFilePathDict.get(blindID+"-02"), 
																							self.trioValidationDict.get( blindID ))
