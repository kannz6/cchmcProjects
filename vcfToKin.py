#!/usr/bin/python
import os
import re
import fileinput

# database connection fields

_n = 0
_fileNumber = 0

def command_line_query(query):
    return "export PGPASSWORD=%s && psql -h %s -d %s -U %s -p %s -c \"%s\"" % (PG_PASSWORD,PG_HOST,PG_DBNAME,PG_USER,PG_PORT,query)

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

		vcfBatchJobsScript = open( "runBatchJobsScript.awk", "w+" )

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

				# vcfBatchJobsScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + batchFileName + "\")}' & \n");
				vcfBatchJobsScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/harvard-vcf/" + batchFileName + "\")}' & \n");

		vcfBatchJobsScript.close()

	def createNPlinkAndKingBatchFiles(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("plinkScriptFilesUsed.txt", "w+")
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(diff-filtered-yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		listOfVcfFileNames = filter( (lambda x : re.match( r'(diff-filtered-harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		batchScript = open("batchPlinkScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			# vcfFile_RegEx = re.search( r'.*(batch_([0-9]+))\.vcf\.gz', vcfFile );#yale
			vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#harvard
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				plinkFileName = "batchPlink_"+ vcfFile_RegEx.group(3) + ".sh"
				with open(plinkFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);
					# batchFileWriter.write("rm -Rf " + vcfFile_RegEx.group(1) + "_kin0\n")#yale
					# batchFileWriter.write("rm -Rf " + vcfFile_RegEx.group(1) + "_kin\n")
					batchFileWriter.write("rm -Rf " + vcfFile_RegEx.group(2) + "_kin\n")#harvard
					# batchFileWriter.write("mkdir " + vcfFile_RegEx.group(1) + "_kin0\n")
					# batchFileWriter.write("mkdir " + vcfFile_RegEx.group(1) + "_kin\n\n")
					batchFileWriter.write("mkdir " + vcfFile_RegEx.group(2) + "_kin\n\n")
					# batchFileWriter.write("rm " + vcfFile_RegEx.group(1) + "/*\n\n")

					# batchFileWriter.write("module load bcftools/1.4\n")
					# filteredVcfFile = "{0}/filtered-{1}".format(vcfFile_RegEx.group(1),vcfFile)#yale
					# filteredVcfFile = "{0}/{1}".format(vcfFile_RegEx.group(2),vcfFile)#harvard
					# batchFileWriter.write("bcftools view -m2 -M2 -v snps {0}/{1} > {2}\n".format(vcfFile_RegEx.group(1), vcfFile, filteredVcfFile))

					batchFileWriter.write("module load plink/1.90b\n")

					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink\n")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+ filteredVcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink\n")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+filteredVcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink")#yale file
					batchFileWriter.write("plink --allow-extra-chr --vcf "+filteredVcfFile+" --make-bed --out "+vcfFile_RegEx.group(2)+"_kin/"+vcfFile_RegEx.group(2)+"_plink")#harvard file
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink")
					# batchFileWriter.write("\ncp {0}_kin/{0}_plink.fam {0}_kin/{1}\n".format(vcfFile_RegEx.group(1),"original_plink.txt"))
					# batchFileWriter.write("\ncp {0}_kin/{0}_plink.fam {0}_kin/{1}\n".format(vcfFile_RegEx.group(2),"original_plink.txt"))

					# batchFileWriter.write("./update_fam {0}\n".format(vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0} {1}".format("update_fam.py", vcfFile_RegEx.group(1)))
					# tmp_fam_file = "{0}/tmp_fam.tx".format()
					# batchFileWriter.write("./update_fam.py {0} > {1}\n".format(vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.fam",vcfFile_RegEx.group(1)+"_kin0/tmp_fam.txt"))
					# batchFileWriter.write("./update_fam.py {0} > {1}\n".format(vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.fam",vcfFile_RegEx.group(1)+"_kin/tmp_fam.txt"))
					# batchFileWriter.write("./update_fam.py {0} > {1}\n".format(vcfFile_RegEx.group(2)+"_kin/"+vcfFile_RegEx.group(2)+"_plink.fam",vcfFile_RegEx.group(2)+"_kin/tmp_fam.txt"))
					# batchFileWriter.write("cp tmp_fam.txt {0}\n".format(vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0}tmp_fam.txt {1}\n".format(vcfFile_RegEx.group(1)+"_kin0/",vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0}tmp_fam.txt {1}\n".format(vcfFile_RegEx.group(1)+"_kin/",vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0}tmp_fam.txt {1}\n".format(vcfFile_RegEx.group(2)+"_kin/",vcfFile_RegEx.group(2)+"_kin/"+vcfFile_RegEx.group(2)+"_plink.fam"))
					batchFileWriter.write("\n\n")
					batchFileWriter.write("module load king/1.4\n")
					# batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.bed --prefix "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+".king\n")
					# batchFileWriter.write("awk '{while(getline){if( $1 !~ /-(01|02)/){print $1 \" \" $1 \" \" $1\"-01 \" $1\"-02 0 -9\\n\" $1 \" \" $1\"-01 0 0 2 -9\\n\" $1 \" \" $1\"-02 0 0 1 -9\"}}}' " + vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.fam > " + vcfFile_RegEx.group(1)+"_kin/"+"temp.txt\n\n")
					# batchFileWriter.write("awk '{i=0;vals[0]=$1;while(getline){ if( $1 !~ /1-.*-(01|02)$/){vals[++i]=$1} };{ for(n=0;n <= i; n++){print vals[n] \" \" vals[n] \" \" vals[n]\"-01 \" vals[n]\"-02 0 -9\\n\" vals[n] \" \" vals[n]\"-01 0 0 2 -9\\n\" vals[n] \" \" vals[n]\"-02 0 0 1 -9\"} } }' " + vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.fam > " + vcfFile_RegEx.group(1)+"/"+"temp.txt\n\n")
					# batchFileWriter.write("cp "+ vcfFile_RegEx.group(1)+"/temp.txt "+ vcfFile_RegEx.group(1)+"/"+ vcfFile_RegEx.group(1)+"_plink.fam\n\n")
					# batchFileWriter.write("cp "+ vcfFile_RegEx.group(1)+"_kin/temp.txt "+ vcfFile_RegEx.group(1)+"_kin/"+ vcfFile_RegEx.group(1)+"_plink.fam\n\n")
					# batchFileWriter.write("module load king/2.0\n")
					# batchFileWriter.write("module load king/1.4\n")
					# batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.bed --prefix "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+".king\n")
					# batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.bed --kinship --prefix "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+".king")
					batchFileWriter.write("king -b "+vcfFile_RegEx.group(2)+"_kin/"+vcfFile_RegEx.group(2)+"_plink.bed --kinship --prefix "+vcfFile_RegEx.group(2)+"_kin/"+vcfFile_RegEx.group(2)+".king")
					# batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.bed --kinship --prefix "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+".king")
					batchFileWriter.write("\n\n")
					# batchFileWriter.write("mv "+vcfFile_RegEx.group(1)+ "* "+vcfFile_RegEx.group(1)+"/")
					# batchFileWriter.write("\n\n")
					batchFileWriter.close()

				# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + plinkFileName + "\")}' & \n");
				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/harvard-vcf/" + plinkFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()


	def diffOrigialVcfAndFilteredVcf(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(.vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("diffScriptFilesUsed.txt", "w+")
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		batchScript = open("batchDiffScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			# vcfFile_RegEx = re.search( r'.*(batch_([0-9]+))\.vcf\.gz', vcfFile );
			# vcfFile_RegEx = re.search( r'.*(yale-(batch_([0-9]+)))\.vcf\.gz', vcfFile );
			vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				diffFileName = "batchDiff_"+ vcfFile_RegEx.group(3) + ".sh"
				with open(diffFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);

					batchFileWriter.write("module load bcftools/1.4\n")
					originalIdFile = "{0}/{0}_filtered_vcf_ids_original.txt".format(vcfFile_RegEx.group(2))
					updatedIdFile = "{0}/{0}_updated_ids.txt".format(vcfFile_RegEx.group(2))
					diffTxtFileName = "{0}/diff_result_{0}.txt" .format(vcfFile_RegEx.group(2))
					batchFileWriter.write("bcftools query -l {0}/{1} > {2}\n\n".format(vcfFile_RegEx.group(2),vcfFile,originalIdFile))
					# command1 = "awk '{print $2}'"
					# command2 = "awk '{if($0 ~ \" \"){print $2} }'"
					# TODO Fix created ids file command 4-20-17
					createdIdsFile = "{0}/created_ids.txt".format(vcfFile_RegEx.group(2))
					createdIdsCommand = "awk '{if ($1 !~ /.*-01$|.*-02$/){print $1 \"\\n\"$1\"-01\\n\"$1\"-02\"} }' {0} > {1}".format(originalIdFile,createdIdsFile)
					batchFileWriter.write("{0}\n".format(createdIdsCommand))
					# TODO Fix skip ids file command 4-20-17
					skipIdsFile = "{0}/{0}_skip_ids.txt".format(vcfFile_RegEx.group(2))
					skipIdsFileCommand = "diff {0} {1} | awk '{if($0 ~ \" \"){print $2} }' > {2}".format(createdIdsFile,originalIdFile,skipIdsFile)
					batchFileWriter.write("{0}\n\n".format(skipIdsFileCommand))
					# batchFileWriter.write("{0} {1}_kin/{1}_plink.fam > {2}\n\n".format(command1,vcfFile_RegEx.group(2),updatedIdFile))
					# batchFileWriter.write("diff {0} {1} | {3} > {2}\n\n".format(originalIdFile,updatedIdFile,diffTxtFileName, command2))
					# diffAndFilteredVcfFile = "{0}/diff-{1}".format(vcfFile_RegEx.group(2),vcfFile)
					# command3 = "bcftools view -m2 -M2 -v snps -S ^"
					# batchFileWriter.write("{3}{0} {4}/{1} > {2}\n".format(diffTxtFileName, vcfFile, diffAndFilteredVcfFile,command3,vcfFile_RegEx.group(2)))
					batchFileWriter.write("\n\n")
					batchFileWriter.close()

				# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + diffFileName + "\")}' & \n");
				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/harvard-vcf/" + diffFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()

	def createRemoveInvalidIdsFile(self):
		listOfRootBatchDirectories = filter( (lambda x : re.match( r'.*(batch_[0-9]+)$', x) ), self.batchDirectories )
		batchScript = open("batchRemoveIdsScript.awk", "w+")
		directoriesUsedChecker = open("createRemoveInvalidIdsFile-dirsUsed.txt", "w+")
		for d in listOfRootBatchDirectories:
			directoriesUsedChecker.write("{0}\n".format(d))
			originalIdFile = "{0}/{0}_filtered_vcf_ids_original.txt".format(d)
			with open("{0}".format(originalIdFile), "r") as originalIdReader:
				originalIdFileContent = originalIdReader.readlines()
			originalIdReader.close()
			originalIdFileContent = [line.strip() for line in originalIdFileContent]

			skipIdsFile = "{0}/{0}_skip_ids.txt".format(d)
			with open("{0}".format(skipIdsFile), "r") as skipIdReader:
				skipIdContent = skipIdReader.readlines()
			skipIdReader.close()
			skipIdContent =  [line.strip() for line in skipIdContent]
			# print fileContent
			# print skipIdContent

			removeIdsFile = "{0}/{0}-remove-ids.txt".format(d)
			removeIdWriter = open("{0}".format(removeIdsFile), "w+")

			for x in originalIdFileContent:
				if originalIdFileContent.count(x[0:7]) == 1 and originalIdFileContent.count(x[0:7]+ "-01") == 1 and originalIdFileContent.count(x[0:7] + "-02") == 1:
					continue
				else:
					# print "{0}".format(x)
					removeIdWriter.write("{0}\n".format(x))
			removeIdWriter.close()
		directoriesUsedChecker.close()

	def reFilterVcfFile(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(.vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("reFilterVcfScriptFilesUsed.txt", "w+")
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		listOfVcfFileNames = filter( (lambda x : re.match( r'(filtered-harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		batchScript = open("batchReFilterScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			# vcfFile_RegEx = re.search( r'.*(batch_([0-9]+))\.vcf\.gz', vcfFile );
			# vcfFile_RegEx = re.search( r'.*(yale-(batch_([0-9]+)))\.vcf\.gz', vcfFile );
			vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				reFilterFileName = "batchReFilter_"+ vcfFile_RegEx.group(3) + ".sh"
				with open(reFilterFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);

					batchFileWriter.write("module load bcftools/1.4\n")
					removeIdsFile = "{0}/{0}-remove-ids.txt".format(vcfFile_RegEx.group(2))	
					diffAndFilteredVcfFile = "{0}/diff-{1}".format(vcfFile_RegEx.group(2),vcfFile)
					bcftoolsFilterCommand = "bcftools view -m2 -M2 -v snps -S ^{0} -Oz".format(removeIdsFile)
					batchFileWriter.write("{0} {1}/{2} > {3}\n".format(bcftoolsFilterCommand, vcfFile_RegEx.group(2), vcfFile, diffAndFilteredVcfFile))
					batchFileWriter.write("\n\n")
					batchFileWriter.close()

				# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + reFilterfFileName + "\")}' & \n");#yale
				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/harvard-vcf/" + reFilterfFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()

	def filterInitialVcfFiles(self,directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("initialVcfFilterScriptFilesUsed.txt", "w+")
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )#yale
		listOfVcfFileNames = filter( (lambda x : re.match( r'(yale-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		# listOfVcfFileNames = filter( (lambda x : re.match( r'(harvard-batch_[0-9]+.vcf.gz)', x) ), listOfVcfFileNames )
		batchScript = open("batchInitialVcfFilterScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			# vcfFile_RegEx = re.search( r'.*(batch_([0-9]+))\.vcf\.gz', vcfFile );#yale
			vcfFile_RegEx = re.search( r'.*(yale-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#harvard
			# vcfFile_RegEx = re.search( r'.*(harvard-(batch_([0-9]+)))\.vcf\.gz', vcfFile );#harvard
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				initialVcfFilterFileName = "batchInitialVcfFilter_"+ vcfFile_RegEx.group(3) + ".sh"
				with open(initialVcfFilterFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);

					batchFileWriter.write("module load bcftools/1.4\n")
					# filteredVcfFile = "{0}/filtered-{1}".format(vcfFile_RegEx.group(1),vcfFile)
					filteredVcfFile = "{0}/filtered-{1}".format(vcfFile_RegEx.group(2),vcfFile)
					# batchFileWriter.write("bcftools view -m2 -M2 -v snps {0}/{1} > {2}\n".format(vcfFile_RegEx.group(1), vcfFile, filteredVcfFile))
					batchFileWriter.write("".format(vcfFile_RegEx.group(2)))
					batchFileWriter.write("bcftools view -m2 -M2 -v snps -Oz {0}/{1} > {2}\n".format(vcfFile_RegEx.group(2), vcfFile, filteredVcfFile))
					batchFileWriter.write("\n\n")
					# batchFileWriter.write("mv {0} {1}".format(vcfFile,vcfFile_RegEx.group(1)+"/"))#yale
					batchFileWriter.write("mv {0} {1}".format(vcfFile,vcfFile_RegEx.group(2)+"/"))#harvard
					batchFileWriter.write("\n\n")
					batchFileWriter.close()

				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + initialVcfFilterFileName + "\")}' & \n");#yale
				# batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/harvard-vcf/" + initialVcfFilterFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()

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
		listOfRootBatchDirectories = filter( (lambda x : re.match( r'.*(batch_[0-9]+)$', x) ), self.batchDirectories )
		# batchScript = open("batchCreateBulkPSQLScript.awk", "w+")
		directoriesUsedChecker = open("CreateBulkPSQLScript.txt", "w+")
		directoriesUsedChecker.write(str(self.batchDirectories))
		for d in listOfRootBatchDirectories:
			directoriesUsedChecker.write("{0}\n".format(d))
			originalIdFile = "{0}/{0}_filtered_vcf_ids_original.txt".format(d)
			outputFile = "{0}/{0}_query_output_.dat".format(d)
			famFile = "{0}/{0}_plink.fam".format(d)
			with open("{0}".format(originalIdFile), "r") as originalIdReader:
				originalIdFileContent = originalIdReader.readlines()
			originalIdReader.close()
			originalIdFileContent = [line.strip() for line in originalIdFileContent]

			triosSqlFile = "{0}/{0}_trios.sql".format(d)
			triosSqlWriter = open("{0}".format(triosSqlFile), "w+")
			triosSqlWriter.write(sqlStatement)
			probands = set()
			for x in originalIdFileContent:
				if originalIdFileContent.count(x[0:7]) != 1 and originalIdFileContent.count(x[0:7]+ "-01") != 1 and originalIdFileContent.count(x[0:7] + "-02") != 1:
					continue
				else:
					if(x != (x[0:7] + "-01") and x != (x[0:7]+ "-02") and x != (x[0:7]+ "-03")):
						probands.add(x)

			for i,p in enumerate(probands):
				if(i == 0):
					if(len(probands) == 1):
						triosSqlWriter.write("'{0}%'\n\t)\n".format(p))
					elif(len(probands) > 1):
						triosSqlWriter.write("'{0}%' or\n".format(p))
				elif(i <= ( len(probands) - 2 )):
					triosSqlWriter.write("\t\tcore_person.blinded_id ILIKE '{0}%' or\n".format(p))
				elif( i == (len(probands) - 1)):
					triosSqlWriter.write("\t\tcore_person.blinded_id ILIKE '{0}%'\n\t)\n".format(p))
			triosSqlWriter.write("ORDER BY core_person.blinded_id;\n")
			triosSqlWriter.close()
			probandWGenderIds = self.getProbandGenders(triosSqlFile, outputFile)
			_ids = map(lambda x: {'id':x[0], 'gender':x[1]}, probandWGenderIds)#get the gender of the proband
			# directoriesUsedChecker.write(str(_ids) + "\n\n____\n\n")
			##########################
			famFileWriter = open(famFile, "w+")
			_usedIds = set()
			for _id in _ids:
			    if ( _id['gender'] == 'M' and _id['id'] not in _usedIds):
					_proband_gender = '1'
					# print "_proband_gender: %s"% _proband_gender
					famFileWriter.write("{0} {0} {0}-01 {0}-02 {1} -9\n{0} {0}-01 0 0 2 -9\n{0} {0}-02 0 0 1 -9\n".format(_id['id'],_proband_gender))
					_usedIds.add(_id['id'])
			    elif( _id['gender'] == 'F' and _id['id'] not in _usedIds):
					_proband_gender = '2'
					# print "_proband_gender: %s"% _proband_gender
					famFileWriter.write("{0} {0} {0}-01 {0}-02 {1} -9\n{0} {0}-01 0 0 2 -9\n{0} {0}-02 0 0 1 -9\n".format(_id['id'],_proband_gender))
					_usedIds.add(_id['id'])
			famFileWriter.close()
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
	    # print os.popen(command).read()
	    # read in and return the query results
	    with open(outputFile,'r') as f:
	        results = [each.strip().split('\t') for each in f.readlines()]

	    return results

	def createBatchJobsFile(self, _inputFile,_vcfFileList, _outDirectory):
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
		_bcftools_command = "module load bcftools/1.4 && bcftools concat -f {0} |  bcftools view -Oz --force-samples -s".format(_vcfFileList)#harvard

		for c,_id in enumerate(inputFileContent):
			global _n; global _fileNumber
			
			if _n == 0:
				batchJobFile.write("{0} {1},{1}-01,{1}-02,".format(_bcftools_command,_id))
				_n = _n + 1
			elif _n < 38 and (_id != inputFileContent[(len(inputFileContent) -1)]) :
				batchJobFile.write("{0},{0}-01,{0}-02,".format(_id))
				_n = _n + 1
			# else _n == 38 or (c == (len(inputFileContent) -1)):
			else:
				# batchJobFile.write("{0},{0}-01,{0}-02 > {1}/batch_{2}/yale-batch_{2}.vcf.gz \n\n".format(_id,_outDirectory,_fileNumber))
				batchJobFile.write("{0},{0}-01,{0}-02 > {1}/batch_{2}/harvard-batch_{2}.vcf.gz \n\n".format(_id,_outDirectory,_fileNumber))
				_n = 0
				_fileNumber = _fileNumber + 1

		batchJobFile.close()