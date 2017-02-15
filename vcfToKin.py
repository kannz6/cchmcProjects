#!/usr/bin/python
import os
import re

class vcfToKin0:

	def __init__(self):
		self.batchFileNames = [];
		self.shellFileBSubCommands = "#!/bin/bash\n#BSUB -W 72:00\n#BSUB -M 8000\n\n"
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
		[ [self.joinPaths(root, f) for f in files] for root, directories, files in os.walk( directory ) ]

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
					awkResultFileWriter.write( awkResultFileLine );
					awkResultFileWriter.write("\n")
				awkResultFileWriter.close()
				if( awkResultFileLine.startswith( "module" ) ):
					self.numberOfFiles = self.numberOfFiles + 1
			# if(n == 0):
				vcfBatchJobsScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + batchFileName + "\")}' & \n");
			# n = n+1

		vcfbatchJobsScript.close()

	def createNPlinkAndKingBatchFiles(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		batchScript = open("batchPlinkScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			vcfFile_RegEx = re.search( r'.*(batch_([0-9]+))\.vcf\.gz', vcfFile );
			if vcfFile_RegEx:
				
				plinkFileName = "batchPlink_"+ vcfFile_RegEx.group(2) + ".sh"
				# print "module load plink/1.90b\nplink --allow-extra-chr --vcf %s --make-bed --out %s"% (vcfFile, vcfFile_RegEx.group(1))
				# if(n == 0):
				with open(plinkFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);
					batchFileWriter.write("mkdir " + vcfFile_RegEx.group(1) + "\n\n")
					batchFileWriter.write("module load plink/1.90b\n")
					batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink")
					batchFileWriter.write("\n\n")
					batchFileWriter.write("module load king/1.4\n")
					batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.bed --prefix "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+".king")
					batchFileWriter.write("\n\n")
					# batchFileWriter.write("mv "+vcfFile_RegEx.group(1)+ "* "+vcfFile_RegEx.group(1)+"/")
					# batchFileWriter.write("\n\n")
					batchFileWriter.close()

				# if(n == 0):
				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + plinkFileName + "\")}' & \n");
					# n = n+1
		batchScript.close()