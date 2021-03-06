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
		self.parentsRelatedCoefficientDict = {}
		self.childRelationCoefficientsDict = {}
		self.sampleGenderCheckDict = {}


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
		# print "ro: %s\nrr: %s"%(ro,rr)
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
		# recordID_regExPart1 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*).*', variablesList[0] );		
		recordID_regExPart1 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*).*', variablesList[1] );
		recordID_regExPart2 = re.search( r'.*\.([0-9]+-[0-9]+?-[0-9]*|[0-9]+-[0-9]+)\..*', variablesList[2] ) 
		# recordID_regExPart3 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[0] ); #king.kin0 format
		recordID_regExPart3 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[1] );#king.kin anf king.kin0 updated fam file format 
		# recordID_regExPart4 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[3] )#king.kin0 and updated fam file format
		recordID_regExPart4 = re.search( r'([A-Z]+-[0-9]+|1-[0-9]+-[0-9]+|1-[0-9]+)', variablesList[2] )#king.kin and updated fam file format


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
			self.appendChildFIDValueList( recordID_regExPart2.group(1) )
		elif recordID_regExPart4:
			recordID_Part4 += recordID_regExPart4.group(1);	bothIDsFound_1_2 = False
			self.appendChildFIDValueList( recordID_regExPart4.group(1) )
		else:
			bothIDsFound_1_2 = False; bothIDsFound_3_4 = False

		# relation_str = variablesList[0] + " " + variablesList[2] + " " + variablesList[7]#king.kin0 format
		# reverse_relation_str = variablesList[2] + " " + variablesList[0] + " " + variablesList[7]#king.kin0 format
		# relation_str = variablesList[1] + " " + variablesList[3] + " " + variablesList[7]#king.kin0 updated fam file format
		# reverse_relation_str = variablesList[3] + " " + variablesList[1] + " " + variablesList[7]#king.kin0 updated fam file format
		# relation_str = variablesList[0] + " " + variablesList[2] + " " + variablesList[8]#king.kin format
		# reverse_relation_str = variablesList[2] + " " + variablesList[0] + " " + variablesList[8]#king.kin format
		relation_str = variablesList[1] + " " + variablesList[2] + " " + variablesList[8]#king.kin format
		reverse_relation_str = variablesList[2] + " " + variablesList[1] + " " + variablesList[8]#king.kin format
		
		self.appendFIDValueList( variablesList[1] ); self.appendFIDValueList( variablesList[2] )#king.kin format append fidvalue
		# self.appendFIDValueList( variablesList[1] ); self.appendFIDValueList( variablesList[3] )#king.kin0 updated fam file format

		# print "Relation Strings: %s\n%s"% (relation_str, reverse_relation_str)
		self.appendRelationStringList( relation_str, reverse_relation_str )
		# only add if the id's have been found so not to throw an error when calling method
		if ( bothIDsFound_1_2 ):
			self.addRelationToDictionary( recordID_Part1, recordID_Part2, relation_str )
		elif ( bothIDsFound_3_4 ):
			self.addRelationToDictionary( recordID_Part3, recordID_Part4, relation_str )
		# else:
		# 	print "Both ids not found!"

	def parseKinFile( self, kinFile ):
		with open( kinFile ) as kf:
			for line in kf:
			#find relevant data in file
				kinFileVariables = line.split()
				if (kinFileVariables[0] == 'FID1' or kinFileVariables[0] == 'FID'):#kin0 and kin format
				# if (kinFileVariables[0] == 'FID'):#kin format
					continue
				# print "******parseKinFile:\n%s"%kinFileVariables
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
		# print "Dict List: %s"% dictList
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

	def setVCFFile( self , fileName ):
		for blinded_id in self.allChildBlindIds:
			if( not self.blindIdFilePathDict.has_key( blinded_id ) ):
				self.blindIdFilePathDict.update( { blinded_id : fileName } )
				self.blindIdFilePathDict.update( { blinded_id + "-01" : fileName } )
				self.blindIdFilePathDict.update( { blinded_id + "-02" : fileName } )
			else:
				print "setVCFFile\n%s-->%s"%(blinded_id, self.blindIdFilePathDict.get(blinded_id))

	def parseDirectory( self ):
		# Run the above function and store its results in a variable.   
		full_file_paths = self.get_filepaths( self.dir )
		# get all the files that match the file extension we are looking for
		# list_of_kin_files = filter( (lambda x : re.match( r'(.*\/.*kin0)', x) ), full_file_paths )
		# list_of_kin_files = filter( (lambda x : re.match( r'(.*\/.*king.kin0)', x) ), full_file_paths )#king.kin0 format
		list_of_kin_files = filter( (lambda x : re.match( r'(.*\/.*king.kin$)', x) ), full_file_paths )#king.kin format#comment out 5-2-17
		# print list_of_kin_files
		# using the list from the kin files, extract the variables we need
		map( self.parseKinFile, list_of_kin_files )#comment out 5-2-17

		# uncomment for input files
		list_of_input_files = filter( (lambda x : re.match( r'(.*\/input_[0-9]+.dat)', x) ), full_file_paths )#comment out 5-2-17
		# print "input files: %s"% list_of_input_files
		map( self.parseExomeFilePaths, list_of_input_files )#comment out 5-2-17

		self.setVCFFile( "/projects/pcgc/prod/data/expedat/yale_vcf/exome_calls.vcf.gz" )
		# self.setVCFFile( "/projects/pcgc/prod/data/expedat/incoming/yale/wide.*.targets.a.f-016.m.vcf.gz" )#comment out 5-2-17

		list_of_genderCheck_files = filter( (lambda x : re.match( r'(.*.sexcheck)', x) ), full_file_paths )
		map( self.setGenderCheckStatuses, list_of_genderCheck_files )
		# self.printGenderDict()

	def parseExomeFilePaths( self, exomeFile ):
		with open( exomeFile, "r" ) as exomeFileReaderO:
			exomeFileContent = exomeFileReaderO.readlines()
			exomeFileReaderO.close()
			exomeFileContent = [line.strip() for line in exomeFileContent] 
			
			for exomeFilePath in exomeFileContent:
				if (not exomeFilePath.startswith("aS")):#uncomment when using input_#.dat files
					continue
				# print exomeFilePath
				# blindIDPath_regEx = re.search( r'.*(1-[0-9]+-[0-9]+|1-[0-9]+).*', exomeFilePath );
				blindIDPath_regEx = re.search( r'.*(1-[0-9]+-[0-9]+|1-[0-9]+).*', exomeFilePath[3:-1] );#use this for parsing input_#.dat file
				if blindIDPath_regEx:
					# print "Blind ID: %s"% blindIDPath_regEx.group(1)
					if ( self.blindIdFilePathDict.has_key( blindIDPath_regEx.group(1) ) ):
						val = self.blindIdFilePathDict.get( blindIDPath_regEx.group(1) )
						# self.blindIdFilePathDict.update( { blindIDPath_regEx.group(1) : val + "," + exomeFilePath } )
						self.blindIdFilePathDict.update( { blindIDPath_regEx.group(1) : val + "," + exomeFilePath[3:-1] } )#use this for parsing input_#.dat file

					else:
						# self.blindIdFilePathDict.update( { blindIDPath_regEx.group(1) : exomeFilePath } )
						self.blindIdFilePathDict.update( { blindIDPath_regEx.group(1) : exomeFilePath[3:-1] } )#use this for parsing input_#.dat file
	
	def validateTrio( self ):
		for r in self.getRelationDictionary():
			idVals = r[0].split() #r[0] contains the FID's/key in the dictionary
			relVals = r[1].split() #r[1] contains the relation/value in the dictionary
			# print "idvals: %s"% idVals
			# print "relvals: %s"% relVals
			print "_______________________________________"
			if ( float(relVals[2] ) > 0.354 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: Duplicate/MZ Twin Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					if( not self.trioValidationDict.has_key( idVals[1] ) ):
						self.trioValidationDict.update( { idVals[1] : "yes" } )
					if( not self.childRelationCoefficientsDict.has_key( idVals[4] ) ):
						self.childRelationCoefficientsDict.update( { idVals[4] : relVals[2] } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[4] : "yes" } )
					if( not self.childRelationCoefficientsDict.has_key( idVals[1] ) ):
						self.childRelationCoefficientsDict.update( { idVals[1] : relVals[2] } )
			elif ( float(relVals[2] ) <= 0.354 and float( relVals[2] ) >= 0.177 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 1st-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					if( not self.trioValidationDict.has_key( idVals[1] ) ):
						self.trioValidationDict.update( { idVals[1] : "yes" } )
					if( not self.childRelationCoefficientsDict.has_key( idVals[4] ) ):
						self.childRelationCoefficientsDict.update( { idVals[4] : relVals[2] } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[4] : "yes" } )
					if( not self.childRelationCoefficientsDict.has_key( idVals[1] ) ):
						self.childRelationCoefficientsDict.update( { idVals[1] : relVals[2] } )
			elif ( float(relVals[2] ) < 0.177 and float( relVals[2] ) >= 0.0884 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 2nd-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[1] : "no" } ); 
					if( not self.childRelationCoefficientsDict.has_key( idVals[4] ) ):
						self.childRelationCoefficientsDict.update( { idVals[4] : relVals[2] } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[4] : "no" } );
					if( not self.childRelationCoefficientsDict.has_key( idVals[1] ) ):
						self.childRelationCoefficientsDict.update( { idVals[1] : relVals[2] } )
			elif ( float(relVals[2] ) < 0.0884 and float( relVals[2] ) >= 0.0442 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 3rd-degree Relation" #+ r[1]
				if ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[1], idVals[4] )
					if( not self.trioValidationDict.has_key( idVals[1] ) ):
						self.trioValidationDict.update( { idVals[1] : "no" } ); 
					
					if( not self.childRelationCoefficientsDict.has_key( idVals[4] ) ):
						self.childRelationCoefficientsDict.update( { idVals[4] : relVals[2] } )
				elif ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ):
					# print "child id: %s, parent id: %s"%( idVals[4], idVals[1] )
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[4] : "no" } ); 

					if( not self.childRelationCoefficientsDict.has_key( idVals[1] ) ):
						self.childRelationCoefficientsDict.update( { idVals[1] : relVals[2] } )
			elif ( float( relVals[2] ) < 0.0442 ):
				print "\n" + r[0] + "\nKinship Value: " + relVals[2] + "\nRelation: 4th-degree Relation" #+ r[1]
				if ( ( ( idVals[1].endswith('-01') or idVals[1].endswith('-02') ) and ( idVals[1][0:7] == idVals[4][0:7] ) ) 
																			and ( ( idVals[4].endswith('-01') or idVals[4].endswith('-02') and ( idVals[4][0:7] == idVals[1][0:7] ) ) ) ):
					if ( ( self.allChildBlindIds.count( idVals[1][0:7] ) == 1 ) and not self.parentsRelatedCoefficientDict.has_key( idVals[1][0:7] ) ):
						self.parentsRelatedCoefficientDict.update( { idVals[1][0:7] : relVals[2] } )
				elif ( ( self.allChildBlindIds.count( idVals[1] ) == 1 and ( idVals[4] == idVals[1] + "-01" or idVals[4] == idVals[1] + "-02" ) ) ):
					if( not self.trioValidationDict.has_key( idVals[1] ) ):
						self.trioValidationDict.update( { idVals[1] : "no" } );

					if( not self.childRelationCoefficientsDict.has_key( idVals[4] ) ):
						self.childRelationCoefficientsDict.update( { idVals[4] : relVals[2] } )
				elif ( ( self.allChildBlindIds.count( idVals[4] ) == 1 and ( idVals[1] == idVals[4] + "-01" or idVals[1] == idVals[4] + "-02" ) ) ):
					if( not self.trioValidationDict.has_key( idVals[4] ) ):
						self.trioValidationDict.update( { idVals[4] : "no" } );

					if( not self.childRelationCoefficientsDict.has_key( idVals[1] ) ):
						self.childRelationCoefficientsDict.update( { idVals[1] : relVals[2] } )

			if ( float(relVals[2] ) > 0.0442 and ( ( ( idVals[1].endswith('-01') or idVals[1].endswith('-02') ) and ( idVals[1][0:7] == idVals[4][0:7] ) )
																			and ( ( idVals[4].endswith('-01') or idVals[4].endswith('-02') ) and ( idVals[4][0:7] == idVals[1][0:7] ) ) ) ):
				if ( ( self.allChildBlindIds.count( idVals[1][0:7] ) == 1 ) and not self.parentsRelatedCoefficientDict.has_key( idVals[1][0:7] ) ):
					self.parentsRelatedCoefficientDict.update( { idVals[1][:-3] : relVals[2] } )

	def printBlindIdRelationData( self ):
		for blindID in self.allChildBlindIds:
			print """Child-id: %s\nChild file: %s\nMom-id: %s | Mom relation value: %s\nMom file: %s\nDad-id: %s | Dad relation value: %s\nDad-file %s \nTrio Validated: %s | Parent-Relation Coefficient: %s\n_____________________________________________________________________________"""%(blindID, self.blindIdFilePathDict.get(blindID), 
			blindID + "-01", self.childRelationCoefficientsDict.get(blindID + "-01"), self.blindIdFilePathDict.get(blindID+"-01"),
			blindID + "-02", self.childRelationCoefficientsDict.get(blindID + "-02"), self.blindIdFilePathDict.get(blindID+"-02"), 
			self.trioValidationDict.get( blindID ),
			self.parentsRelatedCoefficientDict.get( blindID ))

	def getSpreadSheetFortmatedBlindIdRelationData( self ):
		spreadSheetString = "sep=;\nChild_id;Child_file;Mom_id;Mom_relation_value;Mom_file;Dad_id;Dad_relation_value;Dad_file;Trio_Validated;Parent_Relation_Coefficient;\n"
		for blindID in self.allChildBlindIds:
			# global spreadSheetString
			spreadSheetString += blindID + ";" + str( self.blindIdFilePathDict.get( blindID ) ) + ";"
			spreadSheetString += blindID + "-01" + ";" + str( self.childRelationCoefficientsDict.get(blindID + "-01") ) + ";" + str( self.blindIdFilePathDict.get(blindID+"-01") ) + ";"
			spreadSheetString += blindID + "-02" + ";" + str( self.childRelationCoefficientsDict.get(blindID + "-02") ) + ";" + str( self.blindIdFilePathDict.get(blindID+"-02") ) + ";"
			spreadSheetString += str( self.trioValidationDict.get( blindID ) ) + ";" + str( self.parentsRelatedCoefficientDict.get( blindID ) ) + ";;;;;\n"

		return spreadSheetString

	def getSpreadSheetFortmatedBlindIdRelationDataWGender( self ):
		spreadSheetString = "sep=;\nChild_id;Child_gender_status;Child_file;Mom_id;Mom_gender_status;Mom_relation_value;Mom_file;Dad_id;Dad_gender_status;Dad_relation_value;Dad_file;Trio_Validated;Parent_Relation_Coefficient;\n"
		for blindID in self.allChildBlindIds:
			# global spreadSheetString
			spreadSheetString += blindID + ";" + self.sampleGenderCheckDict.get( blindID ) + ";" + str( self.blindIdFilePathDict.get( blindID ) ) + ";"
			spreadSheetString += blindID + "-01" + ";" + self.sampleGenderCheckDict.get( blindID +"-01" ) + ";" + str( self.childRelationCoefficientsDict.get(blindID + "-01") ) + ";" + str( self.blindIdFilePathDict.get(blindID+"-01") ) + ";"
			spreadSheetString += blindID + "-02" + ";" + self.sampleGenderCheckDict.get( blindID +"-02" ) + ";" + str( self.childRelationCoefficientsDict.get(blindID + "-02") ) + ";" + str( self.blindIdFilePathDict.get(blindID+"-02") ) + ";"
			spreadSheetString += str( self.trioValidationDict.get( blindID ) ) + ";" + str( self.parentsRelatedCoefficientDict.get( blindID ) ) + ";;;;;\n"

		return spreadSheetString

	def setGenderCheckStatuses( self, fileName ):
		with open(fileName, "r") as genderCheckFileReader:
			genderCheckFileContent = genderCheckFileReader.readlines()
		genderCheckFileReader.close()
		
		genderCheckFileContent = [line.split() for line in genderCheckFileContent]
		# children = [filter(lambda x: not x.endswith('-01') and not x.endswith('-02'), fileContent ) for fileContent in genderCheckFileContent]
		children = filter(lambda x: '-01' not in x[1] and '-02' not in x[1] , genderCheckFileContent )
		moms = filter(lambda x: '-01' in x[1] , genderCheckFileContent )
		dads = filter(lambda x: '-02' in x[1] , genderCheckFileContent )
		# print "children: {0}\n\nmoms: {1}\n\ndads: {2}\n\n".format(children,moms,dads)

		map(self.setGender,children); map(self.setGender,moms);	map(self.setGender,dads)

	def setGender( self, sampleID ):
		gender = ["Ambiguous", "M", "F" ]
		status = ["SNP-Ambiguous","OK","NOT OK"]
		if(sampleID[0] == 'FID'):
			return

		if self.allChildBlindIds.count( sampleID[1][0:7] ) == 1:
			if not self.sampleGenderCheckDict.has_key( sampleID[1] ):
				if sampleID[4] == "PROBLEM":
					if sampleID[3] == "0":
						self.sampleGenderCheckDict.update( { sampleID[1] : status[0] } )
					else:
						self.sampleGenderCheckDict.update( { sampleID[1] : status[2] } )
				else:
					self.sampleGenderCheckDict.update( { sampleID[1] : "[" + gender[int(sampleID[2])] + "]: " + status[1] } )

	def printGenderDict( self ):
		for k,v in self.sampleGenderCheckDict.items():
			print "[{0}]: {1}".format(k,v)

	def writeBlindedIdsToDoneFile( self, **kwargs ):
		self.allChildBlindIds.sort()
		if kwargs:
			if "o" in kwargs.keys():
				with open("{0}/done.txt".format(kwargs['o']), "w+") as validatedTriosTracker:
					[ validatedTriosTracker.write("{0}\n".format(_id[0:7])) for _id in self.allChildBlindIds ]
		else:
			with open("done.txt", "w+") as validatedTriosTracker:
				[ validatedTriosTracker.write("{0}\n".format(_id[0:7])) for _id in self.allChildBlindIds ]

