#!/usr/bin/python
import fileinput
import math
import os
import random
import re
import stat
import sys


# database connection fields

_n = 0
_fileNumber = 0

#################################
#Steps
#1. createbatchjobsfile
#2. createNVcfBatchFiles
#3. filterInitialVcfFiles
#4. createBulkPSQLScripts
#5. createNPlinkAndKingBatchFiles
#################################

def command_line_query(query):
    return "export PGPASSWORD=%s && psql -h %s -d %s -U %s -p %s -c \"%s\"" % (PG_PASSWORD,PG_HOST,PG_DBNAME,PG_USER,PG_PORT,query)

# http://stackoverflow.com/questions/12791997/how-do-you-do-a-simple-chmod-x-from-within-python
def chmodOfScript(_fileName):#helper function
	_st = os.stat(_fileName); os.chmod(_fileName, _st.st_mode | stat.S_IEXEC)

class vcfToKin0:

	def __init__(self):
		self.batchFileNames = [];
		self.batchDirectories = [];
		self.shellFileBSubCommands = "#!/bin/bash\n#BSUB -W 24:00\n#BSUB -M 10000\n\n"
		self.numberOfFiles = 0

	def joinPaths ( self, rt, fname ):
		filepath = os.path.join( rt, fname )
		self.batchFileNames.append( fname )  # Add it to the list.

	def getFileNames( self, directory ):
		"""
		This function will generate the file names in a directory 
		tree by walking the tree either top-down or bottom-up. For each 
		directory in the tree rooted at directory top (including top itself), 
		it yields a 3-tuple (dirpath, dirnames, filenames).
		"""
		# Walk the tree.
		[ ([self.joinPaths(root, f) for f in files],[self.batchDirectories.append(d) for d in directories]) for root, directories, files in os.walk( directory ) ]

		return self.batchFileNames  # Self-explanatory.

	def createNVcfBatchFiles( self, awkFile ):
		# with open("awkResult.txt", "r") as awkResultFileReader:
		with open(awkFile, "r") as awkResultFileReader:
			awkResultFileContent = awkResultFileReader.readlines()
		awkResultFileReader.close()
		n = 0;
		
		awkResultFileContent = [line.strip() for line in awkResultFileContent]

		vcfBatchJobsScript = open( "runBatchJobsScript.sh", "w+" )

		for awkResultFileLine in awkResultFileContent:
			batchFileName = "batch_"+ str( self.numberOfFiles ) + ".sh"
			# print "file name: %s\ncontent: %s"% (fileName, awkResultFileLine)
			if (awkResultFileLine != ""):
				with open( batchFileName, "w+" ) as awkResultFileWriter:
					awkResultFileWriter.write( self.shellFileBSubCommands );
					awkResultFileWriter.write( "if [ ! -d batch_{0} ]; then mkdir batch_{0}; fi\n".format(str(self.numberOfFiles)) );
					awkResultFileWriter.write( awkResultFileLine );
					awkResultFileWriter.write("\n")
				awkResultFileWriter.close()
				if( awkResultFileLine.startswith( "module" ) ):
					self.numberOfFiles = self.numberOfFiles + 1

				# vcfBatchJobsScript.write("bsub < /scratch/kannz6/temp/vcfFiles/" + batchFileName + " & \n");
				vcfBatchJobsScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + batchFileName + " & \n");

		vcfBatchJobsScript.close()

	def createNPlinkAndKingBatchFiles(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("plinkKingScriptFilesUsed.txt", "w+")
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		batchScript = open("batchPlinkAndKingScript.sh", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			# vcfFile_RegEx = re.search( r'.*(yale-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#yale
			vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#harvard
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				plinkKingFileName = "batchPlinkKing_"+ vcfFile_RegEx.group(3) + ".sh"
				_outDirectory = "{0}/{1}".format(vcfFile_RegEx.group(2),vcfFile_RegEx.group(1) + "_kin")
				_plinkFilesPostFix =  "{0}/{1}".format(_outDirectory,vcfFile_RegEx.group(1) + "_plink")
				_plinkBedFile =  "{0}".format(_plinkFilesPostFix + ".bed")
				_plinkFamFile =  "{0}".format(_plinkFilesPostFix + ".fam")
				_kingFilesPostFix = "{0}/{1}".format(_outDirectory,vcfFile_RegEx.group(1)+".king")
				_tmpFamFile = "tmpFam.txt"
				_gendersTextFile = "{0}/{0}_genders.txt".format(vcfFile_RegEx.group(2),vcfFile_RegEx.group(2))
				_update_fam_program = "update_fam.py"
				with open(plinkKingFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);

					batchFileWriter.write("if [ ! -d {0} ]; then mkdir {0}; else rm {0}/*; fi\n".format(_outDirectory));

					batchFileWriter.write("module load plink/1.90b\n")

					batchFileWriter.write("plink --allow-extra-chr --vcf {0}/{1} --update-sex {2} --check-sex --make-bed --out {3}\n".format(vcfFile_RegEx.group(2),vcfFile,_gendersTextFile,_plinkFilesPostFix))
					batchFileWriter.write("cp {0} {1}\n".format(_update_fam_program, _outDirectory))
					batchFileWriter.write("cd {0}\n".format(_outDirectory))
					batchFileWriter.write("./{0} {1}_plink.fam > {2}\n".format(_update_fam_program,vcfFile_RegEx.group(1),_tmpFamFile))
					batchFileWriter.write("cp {0} {1}_plink.fam\n".format(_tmpFamFile,vcfFile_RegEx.group(1)))
					batchFileWriter.write("cd ../..\n")
					batchFileWriter.write("module load king/1.4\n")
	
					batchFileWriter.write("king -b {0} --kinship --prefix {1}".format(_plinkBedFile,_kingFilesPostFix))
					batchFileWriter.write("\n\n")

					batchFileWriter.close()

				# batchScript.write("bsub < /scratch/kannz6/temp/vcfFiles/" + plinkKingFileName + " & \n");
				batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + plinkKingFileName + " & \n");

		filesUsedChecker.close()
		batchScript.close()

	def filterInitialVcfFiles(self,**kwargs):

		if kwargs:

			directoryFileNames = self.getFileNames( kwargs['directory'] )
			# get all the files that match the file extension we are looking for
			listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
			filesUsedChecker = open("initialVcfFilterScriptFilesUsed.txt", "w+")
			# listOfVcfFileNames = filter( (lambda x : re.match( r'(yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
			listOfVcfFileNames = filter( (lambda x : re.match( r'(harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
			_batchScriptName = "batchInitialVcfFilterScript.sh"
			batchScript = open(_batchScriptName, "w+")
			_filterInitialStorageDir = "filter_initial_vcf_scripts"
			batchScript.write("if [ ! -d {0} ]; then mkdir {0}; else rm {0}/*; fi\nmv batchInitialVcfFilter_* {0}\n".format(_filterInitialStorageDir))
			n = 0
			for i,vcfFile in enumerate(listOfVcfFileNames):
				# vcfFile_RegEx = re.search( r'.*(yale-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#yale
				vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#harvard

				if vcfFile_RegEx:
					filesUsedChecker.write( vcfFile + "\n")
					# initialVcfFilterScriptName = "{0}/{1}.sh".format(kwargs['directory'],"{0}_{1}".format("batchInitialVcfFilter",vcfFile_RegEx.group(3)))
					initialVcfFilterScriptName = "{0}_{1}.sh".format("batchInitialVcfFilter",vcfFile_RegEx.group(3))
					filteredVcfFile = "{0}/filtered-{1}".format(vcfFile_RegEx.group(2),vcfFile)
					_originalIdFile = "{0}/{0}_sample_ids_original.txt".format(vcfFile_RegEx.group(2),vcfFile_RegEx.group(2))

					with open(initialVcfFilterScriptName, "w+") as batchFileWriter:
						batchFileWriter.write(self.shellFileBSubCommands);
						batchFileWriter.write("module load bcftools/1.4\nbcftools query -l {0}/{1} > {2}\n\n".format(vcfFile_RegEx.group(2),vcfFile,_originalIdFile))
						batchFileWriter.write("bcftools view -m2 -M2 -S {0} -v snps -Oz {1}/{2} > {3}\n\n".format(_originalIdFile,vcfFile_RegEx.group(2), vcfFile, filteredVcfFile))
						batchFileWriter.close()

					# batchScript.write("bsub < /scratch/kannz6/temp/vcfFiles/" + initialVcfFilterScriptName + " & \n");#yale
					# batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + initialVcfFilterScriptName + " & \n");
					batchScript.write("bsub < {0}/{1} & \n".format(_filterInitialStorageDir,initialVcfFilterScriptName))

			filesUsedChecker.close()
			batchScript.close()
			chmodOfScript(_batchScriptName)
			# vals = {"directory":/scratch/kannz6/temp/harvard-vcf/"}

	def createBulkPSQLScripts(self, directory):
		sqlStatement = """SELECT core_person.blinded_id AS blinded_id, core_subject.gender AS gender
		FROM pcgc.bulkup_uploadeddatafile
		JOIN pcgc.core_sample
		\tON bulkup_uploadeddatafile.sample_id = core_sample.sample_id
		JOIN pcgc.core_person
		\tON core_sample.person_id = core_person.id
	   	LEFT OUTER JOIN pcgc.core_subject
   		\tON core_person.blinded_id = core_subject.blinded_id
		WHERE
		\t(bulkup_uploadeddatafile.location ILIKE '%p1.fastq.gz' OR bulkup_uploadeddatafile.location ILIKE '%p2.fastq.gz')
		\tAND bulkup_uploadeddatafile.location LIKE '%_EX_%'
		\tAND bulkup_uploadeddatafile.quarantined = false
		\tAND core_person.blinded_id !~ '-.*-0[1-9]'
		\tAND (
		\t\tcore_person.blinded_id ILIKE """
		self.getFileNames( directory )

		listOfRootBatchDirectories = filter( (lambda x : not x.endswith('_kin') and not x.endswith('_kin0') and x.startswith('batch')), self.batchDirectories )
		directoriesUsedChecker = open("createBulkPSQLScriptDirectoriesUsed.txt", "w+")
		directoriesUsedChecker.write(str(listOfRootBatchDirectories) + "\n")
		for d in listOfRootBatchDirectories:
			directoriesUsedChecker.write("{0}\n".format(d))
			originalIdFile = "{0}/{0}_sample_ids_original.txt".format(d)
			outputFile = "{0}/{0}_query_output_.dat".format(d)
			_gendersTextFile = "{0}/{0}_genders.txt".format(d,d)
			with open("{0}".format(originalIdFile), "r") as originalIdReader:
				originalIdFileContent = originalIdReader.readlines()
			originalIdReader.close()
			originalIdFileContent = [line.strip() for line in originalIdFileContent]

			originalIdFileContent = filter(lambda x: "-" in x, originalIdFileContent)
			#from main.py
			children = set(filter(lambda x: not x.endswith('-01') and not x.endswith('-02'), originalIdFileContent))
			moms = set(filter(lambda x: x.endswith('-01'), originalIdFileContent))
			dads = set(filter(lambda x: x.endswith('-02'), originalIdFileContent))
			originalIdFileContent = filter(lambda x: x + '-01' in moms and x + '-02' in dads, children)

			# batchJobFile.write("{0}\n\n".format(str(originalIdFileContent)))
			originalIdFileContent.sort()

			triosSqlFile = "{0}/{0}_trios.sql".format(d)
			triosSqlWriter = open("{0}".format(triosSqlFile), "w+")
			triosSqlWriter.write(sqlStatement)

			for i,p in enumerate(originalIdFileContent):
				if(i == 0):
					if( len(originalIdFileContent) == 1 ):
						triosSqlWriter.write("'{0}%'\n\t)\n".format(p))
					elif( len(originalIdFileContent) > 1 ):
						triosSqlWriter.write("'{0}%' or\n".format(p))
				elif ( i <= (len(originalIdFileContent)-2) ):
					triosSqlWriter.write("\t\tcore_person.blinded_id ILIKE '{0}%' or\n".format(p))
				elif( i == (len(originalIdFileContent)-1) ):
					triosSqlWriter.write("\t\tcore_person.blinded_id ILIKE '{0}%'\n\t)\n".format(p))
			triosSqlWriter.write("ORDER BY core_person.blinded_id;\n")
			triosSqlWriter.close()

			probandWGenderIds = self.getProbandGenders(triosSqlFile, outputFile)
			_ids = map(lambda x: {'id':x[0], 'gender':x[1]}, probandWGenderIds)#get the gender of the proband
			# directoriesUsedChecker.write(str(_ids) + "\n\n____\n\n")
			##########################
			gendersFileWriter = open(_gendersTextFile, "w+")
			_usedIds = set()
			for _id in _ids:
			    if ( _id['gender'] == 'M' and _id['id'] not in _usedIds):
					_proband_gender = '1'
					gendersFileWriter.write("{0} {0} {1}\n{0}-01 {0}-01 2\n{0}-02 {0}-02 1\n".format(_id['id'],_proband_gender))
					_usedIds.add(_id['id'])
			    elif( _id['gender'] == 'F' and _id['id'] not in _usedIds):
					_proband_gender = '2'
					gendersFileWriter.write("{0} {0} {1}\n{0}-01 {0}-01 2\n{0}-02 {0}-02 1\n".format(_id['id'],_proband_gender))
					_usedIds.add(_id['id'])
			gendersFileWriter.close()
			##########################

		directoriesUsedChecker.close()


	def getProbandGenders(self, scriptFile, outputFile):
	    script_file = scriptFile
	    # load the query
	    with open(script_file,'r') as f:
	        query = f.read();
	    # Remove newline characters to permit usage of the client copy command.
	    query = query.replace('\n',' ')
	    # Remove the terminal character because we're going to add code.
	    query = query.replace(';','')
	    # Get rid of comments.
	    query = re.sub('/\*(.*)\*/','',query)
	    # Add the client copy command.
	    query = "\copy (%s) TO \'%s\'" % (query,outputFile)
	    # Convert the query to a call to psql
	    command = command_line_query(query)
	    # run
	    print os.popen(command).read()
	    # read in and return the query results
	    with open(outputFile,'r') as f:
	        results = [each.strip().split('\t') for each in f.readlines()]

	    return results

	def createBatchJobsFile(self, _inputFile,_vcfFileList,_outDirectory):
		batchJobFile = open("{0}/awkResult.txt".format(_outDirectory), "w+")
		with open(_inputFile, "r") as inputFileReader:
			inputFileContent = inputFileReader.readlines()
		inputFileReader.close()
		inputFileContent = [line.strip() for line in inputFileContent]

		inputFileContent = filter(lambda x: "-" in x, inputFileContent)
		#from main.py
		children = set(filter(lambda x: not x.endswith('-01') and not x.endswith('-02'), inputFileContent))
		moms = set(filter(lambda x: x.endswith('-01'), inputFileContent))
		dads = set(filter(lambda x: x.endswith('-02'), inputFileContent))
		inputFileContent = filter(lambda x: x + '-01' in moms and x + '-02' in dads, children)


		# batchJobFile.write("{0}\n\n".format(str(inputFileContent)))
		inputFileContent.sort()

		# _bcftools_command = "module load bcftools/1.4 && bcftools view -Oz {0} --force-samples -s".format(_vcfFileList)#yale
		# _bcftools_command = "module load bcftools/1.4 && bcftools view -Oz {0} -s".format(_vcfFileList)#yale
		# _bcftools_command = "module load bcftools/1.4 && bcftools concat -f {0} |  bcftools view -Oz --force-samples -s".format(_vcfFileList)#harvard
		_bcftools_command = "module load bcftools/1.4 && bcftools concat -f {0} |  bcftools view -Oz -s".format(_vcfFileList)#harvard
		for c,_id in enumerate(inputFileContent):
			global _n; global _fileNumber
			if _n == 0:
				batchJobFile.write("{0} {1},{1}-01,{1}-02,".format(_bcftools_command,_id))
				_n = _n + 1
			elif _n < 38 and (_id != inputFileContent[(len(inputFileContent) -1)]) :
				batchJobFile.write("{0},{0}-01,{0}-02,".format(_id))
				_n = _n + 1
			else:
				# batchJobFile.write("{0},{0}-01,{0}-02 > {1}/batch_{2}/yale-batch_{2}.vcf.gz \n\n".format(_id,_outDirectory,_fileNumber))
				batchJobFile.write("{0},{0}-01,{0}-02 > {1}/batch_{2}/harvard-batch_{2}.vcf.gz \n\n".format(_id,_outDirectory,_fileNumber))
				_n = 0
				_fileNumber = _fileNumber + 1
		batchJobFile.close()

	def createBatchConcordanceFiles(self, **kwargs):
		####-Params
		####---{in=<input text file of vcfs to loop through>, -S=<sample-ids txt file>, f1=<file one>, f2=<file two>, out=<full path of output file to write results to>, 
		# command=<stats|gtcheck>,outDir=<directory to place output files>,g=<reference genotypes to compare against>}
		####---i.e. bcftools stats -S yale-sample-ids.txt ../exome_calls.vcf.gz ../harvard.vcf.gz > concordance-yale-harvard-4-28-17.txt
		_batchScriptName = "batchConcordanceCreationScript.sh"
		batchScript = open(_batchScriptName, "w+")
		_bulkConcordanceScriptsStorageDir = "concordance_scripts"
		batchScript.write("if [ ! -d {0} ]; then mkdir {0}; else rm {0}/*; fi\nmv batchConcordanceCreation_* {0}\n".format(_bulkConcordanceScriptsStorageDir))
		_moduleLoadCommand = "module load "
		_concordanceScript = "{0}".format(self.shellFileBSubCommands)
		_vcfFilePathsDict = {}
		_bulkFilePathsScriptsList = []
		if kwargs:
			if'outDir' in kwargs.keys():
				batchScript.write("if [ ! -d {0} ]; then mkdir {0}; fi\n".format(kwargs['outDir']))
			
			if 'bcfv' in kwargs.keys():
				_moduleLoadCommand += "{0}\n".format(kwargs['bcfv'])
			else:
				_moduleLoadCommand += "{0}\n".format("bcftools/1.4")

			if 'htslib' in kwargs.keys():
				_moduleLoadCommand += "module load {0}\n".format(kwargs['htslib'])
			else:
				_moduleLoadCommand += "module load {0}\n".format("htslib/1.3")

			_bcftool = ""
			if 'command' in kwargs.keys():
				_bcftool = kwargs['command']
			else:
				_bcftool = 'stats'

			if 'xS' in kwargs.keys():
				_sampleIdFile = kwargs['xS']
				_concordanceScript += "{0}\nbcftools {1} -S ^{2} ".format(_moduleLoadCommand,_bcftool,_sampleIdFile)
			elif 'S' in kwargs.keys():
				_sampleIdFile = kwargs['S']
				_concordanceScript += "{0}\nbcftools {1} -S {2} ".format(_moduleLoadCommand,_bcftool,_sampleIdFile)
			else:
				_concordanceScript += "{0}\nbcftools {1} ".format(_moduleLoadCommand,_bcftool)

			if 'stats' in _bcftool:
				if 'e' in kwargs.keys():
					_exclude = kwargs['e']
					_concordanceScript += "-e \'{0}\' ".format(kwargs['e'])
				elif 'i' in kwargs.keys():
					_include = kwargs['i']
					_concordanceScript += "-i \'{0}\' ".format(kwargs['i'])
			elif 'gtcheck' in _bcftool:
				if 'g' in kwargs.keys():
					_genotypes = kwargs['g']
					_concordanceScript += "-g {0} ".format(_genotypes)
			
			_file1 = kwargs['f1']; _concordanceScript += "{0} ".format(_file1)

			if 'in' in kwargs.keys():
				_inputTextFileOfVcfFilePaths = kwargs['in']
				with open(_inputTextFileOfVcfFilePaths, 'r') as _inputTextFileOfVcfFilePathsReader:
					_inputTextFileOfVcfFilePathsContent = _inputTextFileOfVcfFilePathsReader.readlines()
				_inputTextFileOfVcfFilePathsReader.close()
				_inputTextFileOfVcfFilePathsContent = [line.strip() for line in _inputTextFileOfVcfFilePathsContent]
				[ _vcfFilePathsDict.update({ i+1 : v }) for i,v in enumerate(_inputTextFileOfVcfFilePathsContent) ]

				[_bulkFilePathsScriptsList.append("{0} {1} > {2}/{3}_concordance.txt\n\nexit\n\n".format(_concordanceScript,v,kwargs['outDir'],k)) for k,v in _vcfFilePathsDict.items() ]
				
				for c,e in enumerate(_bulkFilePathsScriptsList):
					bulkConcordanceFileName = "batchConcordanceCreation_{0}.sh".format(c+1)

					with open(bulkConcordanceFileName, "w+") as batchFileWriter:
						batchFileWriter.write("{0}".format(e))
					batchFileWriter.close()

					# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + bulkConcordanceFileName + "\")}' & \n");#yale
					batchScript.write("bsub < {0}/{1} & \nsleep 1".format(_bulkConcordanceScriptsStorageDir,bulkConcordanceFileName))
					# batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + bulkConcordanceFileName + " & \nsleep 2\n");#harvard-dir

			elif 'f2' in kwargs.keys() or 'gtcheck' in _bcftool:
				if 'gtcheck' not in _bcftool:
					_file2 = kwargs['f2']; _concordanceScript += "{0} ".format(_file2)
				_concordanceResultsFile = "{0}/{1}".format(kwargs['outDir'],"concordance-result.txt"); _concordanceScript += "> {0}\n\nexit\n\n".format(_concordanceResultsFile)
				bulkConcordanceFileName = "batchConcordanceCreation_1.sh"
				with open(bulkConcordanceFileName, "w+") as batchFileWriter:
					batchFileWriter.write("{0}".format(_concordanceScript))
				batchFileWriter.close()
				# _bulkConcordanceScriptsStorageDir
				batchScript.write("bsub < {0}/{1} & \n".format(_bulkConcordanceScriptsStorageDir,bulkConcordanceFileName))
				# batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + bulkConcordanceFileName + " &\n");#harvard-dir

		batchScript.close()
		chmodOfScript(_batchScriptName)
		# vals= {"S": "4-28-17/yale-sample-ids.txt","outDir" : "5-15-17", "bcfv" : "bcftools/1.3", "in" : "harvard.txt", "f1":"exome_calls.vcf.gz","e":r"FMT/QD<5 | FMT/AB < 0.2 | FMT/GQ < 10 | FMT/DP < 8"}
		# createBatchConcordanceFiles(**vals)

	def createBatchGrepScript(self, **kwargs):
		####-kw arguments to search for in passed in variables dictionary
		####---{inDir=<input directory where files to be greped reside>, fileExt=<extension to include in files to be grepped>,g=<what to grep in file(s)format: {"1":x,"2":y..."n":z}>, outDir=<directory to store output results>, outName=<output file name to append to result(s)> }
		####-create bulk scripts that grep parameter in file
		####---i.e. 
		# get all the files that match the file extension we are looking for
		#i.e. cat 5-15-17/96_concordance.txt | grep "GCsS"
		_grepFileName = "grepScript.sh"
		batchScript = open(_grepFileName, "w+")
		if kwargs:
			directoryFileNames = self.getFileNames( kwargs['inDir'] )
			if "ext" in kwargs.keys():
				_ext = "{0}".format(kwargs['fileExt'])
			else:
				_ext = "concordance-result.txt"
			grepCommand = "{0}\nif [ ! -d {1}/{2} ]; then mkdir {1}/{2}; else rm {1}/{2}/*; fi\ncat {1}/*{3} | egrep \'".format(self.shellFileBSubCommands,kwargs['inDir'],kwargs['outDir'],_ext)
			for i,v in enumerate(kwargs['g']):
				if i == (len(kwargs['g']) - 1) :
					grepCommand += "{0}".format(kwargs['g'][v])
				else:
					grepCommand += "{0}|".format(kwargs['g'][v])

			grepCommand += "\' > {0}\nexit\n\n".format("{0}/{1}/{2}".format(kwargs['inDir'],kwargs['outDir'],kwargs['outName']))
			
			with open(_grepFileName, "w+") as grepFileWriter:
				grepFileWriter.write("{0}".format(grepCommand))
			grepFileWriter.close()

			chmodOfScript(_grepFileName)

		# grepVals = {"inDir" : "5-15-17", "fileExt" : ".txt", "g" : {"1" : "GCsS", "2" : "GCiS"}, "outDir" : "yale-harvard-concordance", "outName" : "snpsAndIndels-yale-harvard.txt"}

	def bcfMergeConcatVcfFiles(self, **kwargs):
		# bcftools concat -Oz -f harvard.txt -o /scratch/kannz6/temp/harvard-vcf/harvard-exome-5-17-17.vcf.gz
		# bcftools index harvard-exome-5-17-17.vcf.gz
		#kwargs = {"out":<filename of output>,"bcfv":<version of bcftools>,"htslibv":<htslib version>,"outDir":<directory to write output file to>}
		if kwargs:
			_listOfVcfFiles = kwargs['inFile']

			if "out" in kwargs.keys():
				if "vcf.gz" in kwargs['out']:
					_outFile = kwargs['out']
				else:
					_outFile = "{0}.vcf.gz".format(kwargs['out'])

			_command = "{0}".format(self.shellFileBSubCommands)
			if "bcfv" in kwargs.keys() and kwargs['bcfv']:
				if "bcftools" in kwargs['bcfv']:
					if "/" in kwargs['bcfv']:
						_command += "module load {0}\n".format(kwargs['bcfv'])
					else:
						_command += "module load /{0}\n".format(kwargs['bcfv'])
				else:
					_command += "module load bcftools/{0}\n".format(kwargs['bcfv'])
			else:
				_command += "module load {0}\n".format("bcftools/1.3")
			
			if "htslibv" in kwargs.keys() and kwargs['htslibv']: 
				if "htslib" in kwargs['htslibv']:
					if "/" in kwargs['htslibv']:
						_command += "module load {0}\n".format(kwargs['htslibv'])
					else:
						_command += "module load /{0}\n".format(kwargs['htslibv'])
				else:
					_command += "module load htslib/{0}\n".format(kwargs['htslibv'])
			else:
				_command += "module load {0}\n".format("htslib/1.3")

			if "outDir" in kwargs.keys() and kwargs['outDir']:
				_command += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format(kwargs['outDir'])
				_outFile = "{0}/{1}".format(kwargs['outDir'],_outFile)

			if "concat" in kwargs['tool']:
				_command += "bcftools {0} -Oz -f {1} -o {2}\n".format(kwargs['tool'],_listOfVcfFiles,_outFile)
			elif "merge" in kwargs['tool']:
				_vcfFilePathsDict = {}
				with open(_listOfVcfFiles, 'r') as _inputTextFileOfVcfFilePathsReader:
					_inputTextFileOfVcfFilePathsContent = _inputTextFileOfVcfFilePathsReader.readlines()
				_inputTextFileOfVcfFilePathsReader.close()
				_inputTextFileOfVcfFilePathsContent = [line.strip() for line in _inputTextFileOfVcfFilePathsContent]
				[ _vcfFilePathsDict.update({ i+1 : v }) for i,v in enumerate(_inputTextFileOfVcfFilePathsContent) ]
				_indexCommands = [ "bcftools index -f {0}\n".format(_vcfFile) for _vcfFile in _vcfFilePathsDict.values() ]
				for ic in _indexCommands:
					_command += ic
				_command += "bcftools {0} -Oz -l {1} --force-samples -o {2}\n".format(kwargs['tool'],_listOfVcfFiles,_outFile)
			_command += "bcftools index {0}\n\nexit\n\n".format(_outFile)

			concatScriptFileName = "mergeConcatVcfsScript.sh"
			with open(concatScriptFileName, "w+") as fileWriter:
				fileWriter.write("{0}".format(_command))
			fileWriter.close()
			# batchScript.write("bsub < /scratch/kannz6/temp/vcfFiles/" + concatScriptFileName &\n");#yale


	def testBulkNoAwk(self, **kwargs):
		####-kw arguments to search for in passed in variables dictionary
		####---{inDir=<input directory where files to be greped reside>, fileExt=<extension to include in files to be grepped>,g=<what to grep in file(s)format: {"1":x,"2":y..."n":z}>, outDir=<directory to store output results>, outName=<output file name to append to result(s)> }
		####-create bulk scripts that grep parameter in file
		####---i.e. 
		# get all the files that match the file extension we are looking for
		#i.e. cat 5-15-17/96_concordance.txt | grep "GCsS"
		_batchScriptName = "testBulkNoAwk.sh"
		batchScript = open(_batchScriptName, "w+")
		_moduleLoadCommand = "module load "
		_concordanceScript = "echo \"hello world, didn't use awk for these jobs!!!\""
		_vcfFilePathsDict = {}
		_bulkFilePathsScriptsList = []
		if kwargs:
			if kwargs['in']:
				_inputTextFileOfVcfFilePaths = kwargs['in']
			
			with open(_inputTextFileOfVcfFilePaths, 'r') as _inputTextFileOfVcfFilePathsReader:
				_inputTextFileOfVcfFilePathsContent = _inputTextFileOfVcfFilePathsReader.readlines()
			_inputTextFileOfVcfFilePathsReader.close()
			_inputTextFileOfVcfFilePathsContent = [line.strip() for line in _inputTextFileOfVcfFilePathsContent]
			[ _vcfFilePathsDict.update({ i+1 : v }) for i,v in enumerate(_inputTextFileOfVcfFilePathsContent) ]
			[_bulkFilePathsScriptsList.append("{0}{1} > {2}/{3}_noAwkTestResult.txt\n\nexit\n\n".format(self.shellFileBSubCommands,_concordanceScript,kwargs['outDir'],k)) for k,v in _vcfFilePathsDict.items() ]
				
			for c,e in enumerate(_bulkFilePathsScriptsList):
				bulkConcordanceFileName = "bulkTestNoAwk_{0}.sh".format(c+1)

				with open(bulkConcordanceFileName, "w+") as batchFileWriter:
					batchFileWriter.write("{0}".format(e))
				batchFileWriter.close()

				# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + bulkConcordanceFileName + "\")}' & \n");#yale
				# batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + bulkConcordanceFileName + " & \nsleep 2\n");#harvard-dir
				if c < (len(_bulkFilePathsScriptsList) - 1):
					batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + bulkConcordanceFileName + " &\nsleep 1\n");#harvard-dir
				else:
					batchScript.write("bsub < /scratch/kannz6/temp/harvard-vcf/" + bulkConcordanceFileName + " &\n");#harvard-dir

		batchScript.close()
		chmodOfScript(_batchScriptName)
		# tvals= {"outDir" : "test", "in" : "harvard.txt"}

	def createVcfFileSampleIdsFile(self, **kwargs):#kwargs: i,d,bcfv
		if kwargs:
			script = open("getVcfFileSampleIdsScript.sh", "w+")
			cmd = "{0}".format(self.shellFileBSubCommands)

			if "i" not in kwargs.keys():
				sys.exit("No input vcf file in dictionary!")
			else:
				if "d" not in kwargs.keys():
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format("temp")
					_d = "temp"
				else:
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format(kwargs['d'])
					_d = kwargs['d']

				_vcfSampleIds = "{0}/{1}".format(_d,"sample_ids.txt")

				if "bcfv" not in kwargs.keys():
					cmd += "module load bcftools/1.4\nbcftools query -l {0} > {1}\n\nexit\n".format(kwargs['i'], _vcfSampleIds)
				else:
					cmd += "module load bcftools/{0}\nbcftools query -l {1} > {2}\n\nexit\n".format(kwargs['bcfv'],kwargs['i'],_vcfSampleIds)
			script.write(cmd)
			script.close()
			# chmodOfScript(script)
			# vals = {"i":"harvard-exome-5-17-17.vcf.gz","d":"5-25-17-refilterVcf"}
			# batchRunner.createVcfFileSampleIdsFile(**vals)
			
	def filterSampleIds(self, **kwargs):#kwargs: s,d
		if kwargs:

			if "s" not in kwargs.keys():
				sys.exit("No sample ids file in dictionary!")
			with open(kwargs['s'], "r") as inputFileReader:
				inputFileContent = inputFileReader.readlines()
			inputFileReader.close()
			inputFileContent = [line.strip() for line in inputFileContent]

			inputFileContent = filter(lambda x: "-" in x, inputFileContent)
			#from main.py
			children = set(filter(lambda x: not x.endswith('-01') and not x.endswith('-02'), inputFileContent))
			moms = set(filter(lambda x: x.endswith('-01'), inputFileContent))
			dads = set(filter(lambda x: x.endswith('-02'), inputFileContent))
			inputFileContent = filter(lambda x: x + '-01' in moms and x + '-02' in dads, children)

			inputFileContent.sort()
			if "d" in kwargs.keys():
				_outputFile = "{0}/filtered-sample-ids.txt".format(kwargs['d'])
			else:
				_outputFile = "{0}/filtered-sample-ids.txt".format("temp")
			updatedSamplesFile = open(_outputFile, "w+")

			if "subset" not in kwargs.keys():
				[ updatedSamplesFile.write("{0}\n{0}-01\n{0}-02\n".format(_id)) for _id in inputFileContent ]
			else:
				ss = []; children = list(children)
				if "using-count" not in kwargs["subset"].keys() and "size" not in kwargs["subset"].keys():
					[ ss.append(children[c]) for c in range(0,len(children),random.randint(1,len(children) - 1)) ]			
				elif "using-count" in kwargs["subset"].keys() and isinstance(kwargs['subset']['using-count'],int):
					[ ss.append(children[c]) for c in range(0,len(children),int(kwargs['subset']['using-count'])) ]
				elif "size" in kwargs["subset"].keys() and isinstance(kwargs['subset']['size'],int) and kwargs['subset']['size'] <= len(inputFileContent) - 1:				
					[ ss.append(children[random.randint(0,len(children)-1)]) for c in range(0,len(children)) if len(ss) * 3 <= kwargs['subset']['size'] and children[random.randint(0,len(children)-1)] not in ss]

				inputFileContent = filter(lambda x: x + '-01' in moms and x + '-02' in dads, ss);

				if "sort" in kwargs['subset'].keys() and kwargs['subset']['sort'] is True:
					inputFileContent.sort()

				[ updatedSamplesFile.write("{0}\n{0}-01\n{0}-02\n".format(_id)) for _id in inputFileContent ]

			updatedSamplesFile.close()
			
			# vals = {"s":"5-25-17-refilterVcf/sample_ids.txt","d":"5-25-17-refilterVcf"}
			# batchRunner.filterSampleIds(**vals)

	def reFilterSingleVCF(self, **kwargs):#kwargs: i,d,bcfv,s,o,htslibv
		if kwargs:
			script = open("reFilterSingleVCFScript.sh", "w+")
			cmd = "{0}".format(self.shellFileBSubCommands)

			if "i" not in kwargs.keys():
				sys.exit("No input vcf file in dictionary!")
			else:
				if "d" not in kwargs.keys():
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format("temp")
					_d = "temp"
				else:
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format(kwargs['d'])
					_d = kwargs['d']

				_vcfSampleIds = "{0}/{1}".format(_d,"sample_ids.txt")

				if "bcfv" not in kwargs.keys():
					cmd += "module load bcftools/1.4\n"
				else:
					cmd += "module load bcftools/{0}\n".format(kwargs['bcfv'])

				if "s" not in kwargs.keys():
					sys.exit("No sample ids file in dictionary!")
				else:
					if "o" in kwargs.keys():
						_outFile = kwargs['o']
						cmd += "bcftools view -m2 -M2 -S {0} -v snps -Oz {1} > {2}/{3}\n".format(kwargs['s'], kwargs['i'], kwargs['d'], kwargs['o'])
					else:
						_outFile = "{0}/filtered-{1}".format(_d,kwargs['i'])
						cmd += "bcftools view -m2 -M2 -S {0} -v snps -Oz {1} > {2}\n".format(kwargs['s'], kwargs['i'], _outFile)

				if "htslibv" not in kwargs.keys():
					cmd += "module load htslib/1.3\nbcftools index {0}\n\nexit\n".format(_outFile)
				else:
					cmd += "module load htslib/{0}\nbcftools index {1}\n\nexit\n".format(kwargs['htslibv'],_outFile)

			script.write(cmd)
			script.close()
			# chmodOfScript(script)
			# vals = {"i":"harvard-exome-5-17-17.vcf.gz",s":"5-25-17-refilterVcf/filtered-sample-ids.txt","d":"5-25-17-refilterVcf"}
			# batchRunner.reFilterSingleVCF(**vals)

	def computeIntersection(self, **kwargs):#vcf-a,vcf-b,o,f-ex,f-in,d,htslibv,bcfv
		if kwargs:
			script = open("computeIntersection.sh", "w+")
			cmd = "{0}".format(self.shellFileBSubCommands)

			if "vcf-a" not in kwargs.keys() or "vcf-b" not in kwargs.keys():
				sys.exit("Missing input vcf file(s). Must supply vcf-a AND vcf-b.")
			else:
				if "d" not in kwargs.keys():
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format("temp")
					_d = "temp"
				else:
					cmd += "if [ ! -d {0} ]; then mkdir {0}; fi\n".format(kwargs['d'])
					_d = kwargs['d']

				if "htslibv" not in kwargs.keys():
					cmd += "module load htslib/1.3\n"
				else:
					cmd += "module load htslib/{0}\n".format(kwargs['htslibv'])

				if "bcfv" not in kwargs.keys():
					cmd += "module load bcftools/1.4\n"
				else:
					cmd += "module load bcftools/{0}\n".format(kwargs['bcfv'])

				cmd += "bcftools isec -p {0} -Oz ".format(kwargs['d'])

				if "f-in" in kwargs.keys():
					cmd += "-i \'{0}\' ".format(kwargs['f-in'])

	 			if "f-ex" in kwargs.keys():
					cmd += "-e \'{0}\' ".format(kwargs['f-ex'])

				cmd += "{0} {1}\n".format(kwargs['vcf-a'], kwargs['vcf-b'])

			script.write(cmd)
			script.close()
		else:
			sys.exit("You must provide an arguments dictionary to computeIntersection:\nargs: vcf-a, vcf-b, o, f-ex, f-in, htslib, bcfv, d")
		# vals = {"d":"6-1-17-intersection",vcf-a":"filtered-harvard-exome-5-17-17.vcf.gz","vcf-b":"exome_calls.vcf.gz",}
		# batchRunner.computeIntersection(**vals)