#!/usr/bin/python
import os
import re

class vcfToKin0:

	def __init__(self):
		self.batchFileNames = [];
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

				vcfBatchJobsScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + batchFileName + "\")}' & \n");

		vcfbatchJobsScript.close()

	def createNPlinkAndKingBatchFiles(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("plinkScriptFilesUsed.txt", "w+")

		batchScript = open("batchPlinkScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			vcfFile_RegEx = re.search( r'.*diff-filtered-(batch_([0-9]+))\.vcf\.gz', vcfFile );
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				plinkFileName = "batchPlink_"+ vcfFile_RegEx.group(2) + ".sh"
				with open(plinkFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);
					batchFileWriter.write("rm -Rf " + vcfFile_RegEx.group(1) + "_kin0\n")
					batchFileWriter.write("rm -Rf " + vcfFile_RegEx.group(1) + "_kin\n")
					batchFileWriter.write("mkdir " + vcfFile_RegEx.group(1) + "_kin0\n")
					batchFileWriter.write("mkdir " + vcfFile_RegEx.group(1) + "_kin\n\n")
					# batchFileWriter.write("rm " + vcfFile_RegEx.group(1) + "/*\n\n")

					# batchFileWriter.write("module load bcftools/1.3\n")
					filteredVcfFile = "{0}/{1}".format(vcfFile_RegEx.group(1),vcfFile)
					# batchFileWriter.write("bcftools view -m2 -M2 -v snps {0} > {1}\n".format(vcfFile, filteredVcfFile))

					batchFileWriter.write("module load plink/1.90b\n")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink\n")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+ filteredVcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink\n")
					batchFileWriter.write("plink --allow-extra-chr --vcf "+filteredVcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink")
					# batchFileWriter.write("plink --allow-extra-chr --vcf "+vcfFile+" --make-bed --out "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink")
					batchFileWriter.write("\n")
					# batchFileWriter.write("./update_fam {0}\n".format(vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0} {1}".format("update_fam.py", vcfFile_RegEx.group(1)))
					# tmp_fam_file = "{0}/tmp_fam.tx".format()
					# batchFileWriter.write("./update_fam.py {0} > {1}\n".format(vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.fam",vcfFile_RegEx.group(1)+"_kin0/tmp_fam.txt"))
					batchFileWriter.write("./update_fam.py {0} > {1}\n".format(vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.fam",vcfFile_RegEx.group(1)+"_kin/tmp_fam.txt"))
					# batchFileWriter.write("cp tmp_fam.txt {0}\n".format(vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					# batchFileWriter.write("cp {0}tmp_fam.txt {1}\n".format(vcfFile_RegEx.group(1)+"_kin0/",vcfFile_RegEx.group(1)+"_kin0/"+vcfFile_RegEx.group(1)+"_plink.fam"))
					batchFileWriter.write("cp {0}tmp_fam.txt {1}\n".format(vcfFile_RegEx.group(1)+"_kin/",vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.fam"))
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
					batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+"_plink.bed --kinship --prefix "+vcfFile_RegEx.group(1)+"_kin/"+vcfFile_RegEx.group(1)+".king")
					# batchFileWriter.write("king -b "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+"_plink.bed --kinship --prefix "+vcfFile_RegEx.group(1)+"/"+vcfFile_RegEx.group(1)+".king")
					batchFileWriter.write("\n\n")
					# batchFileWriter.write("mv "+vcfFile_RegEx.group(1)+ "* "+vcfFile_RegEx.group(1)+"/")
					# batchFileWriter.write("\n\n")
					batchFileWriter.close()

				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + plinkFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()


	def diffVcfAndFamIds(self, directory):

		directoryFileNames = self.getFileNames( directory )
		# get all the files that match the file extension we are looking for
		listOfVcfFileNames = filter( (lambda x : re.match( r'.*(vcf.gz)', x) ), directoryFileNames )
		filesUsedChecker = open("diffScriptFilesUsed.txt", "w+")

		batchScript = open("batchDiffScript.awk", "w+")
		n = 0
		for vcfFile in listOfVcfFileNames:
			vcfFile_RegEx = re.search( r'.*filtered-(batch_([0-9]+))\.vcf\.gz', vcfFile );
			if vcfFile_RegEx:
				filesUsedChecker.write( vcfFile + "\n")
				diffFileName = "batchDiff_"+ vcfFile_RegEx.group(2) + ".sh"
				with open(diffFileName, "w+") as batchFileWriter:
					batchFileWriter.write(self.shellFileBSubCommands);

					batchFileWriter.write("module load bcftools/1.4\n")
					originalIdFile = "{0}/{0}_filtered_vcf_ids_original.txt".format(vcfFile_RegEx.group(1))
					updatedIdFile = "{0}/{0}_updated_ids.txt".format(vcfFile_RegEx.group(1))
					diffTxtFileName = "{0}/diff_result_{0}.txt" .format(vcfFile_RegEx.group(1))
					batchFileWriter.write("bcftools query -l {0} > {1}\n\n".format(vcfFile,originalIdFile))
					command1 = "awk '{print $2}'"
					command2 = "awk '{if($0 ~ " "){print $2} }'"
					batchFileWriter.write("{2} {0}_kin/{0}_plink.fam > {1}\n\n".format(vcfFile_RegEx.group(1),updatedIdFile, command1))
					batchFileWriter.write("diff {0} {1} | {3} > {2}\n\n".format(originalIdFile,updatedIdFile,diffTxtFileName, command2))
					diffAndFilteredVcfFile = "{0}/diff-{1}".format(vcfFile_RegEx.group(1),vcfFile)
					command3 = "bcftools view -m2 -M2 -v snps -S ^"
					batchFileWriter.write("{3}{0} {1} > {2}\n".format(diffTxtFileName, vcfFile, diffAndFilteredVcfFile,command3))
					batchFileWriter.write("\n\n")
					batchFileWriter.close()

				batchScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + diffFileName + "\")}' & \n");

		filesUsedChecker.close()
		batchScript.close()