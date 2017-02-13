#!/usr/bin/python

with open("awkResult.txt", "r") as awkResultFileReader:
	awkResultFileContent = awkResultFileReader.readlines()
awkResultFileReader.close()

n = 0;
shellFileBSubCommands = "#!/bin/bash\n#BSUB -W 72:00\n#BSUB -M 8000\n\n"
awkResultFileContent = [line.strip() for line in awkResultFileContent]

batchJobsScript = open("runBatchJobsScript.awk", "w+")

for awkResultFileLine in awkResultFileContent:
	fileName = "batch_"+ str(n) + ".sh"
	# print "file name: %s\ncontent: %s"% (fileName, awkResultFileLine)
	if (awkResultFileLine != ""):
		with open(fileName, "w+") as awkResultFileWriter:
			awkResultFileWriter.write(shellFileBSubCommands);
			awkResultFileWriter.write(awkResultFileLine);
			awkResultFileWriter.write("\n")
		awkResultFileWriter.close()
		if(awkResultFileLine.startswith("module")):
			n = n+1
	# if(n == 0):
		batchJobsScript.write("awk 'BEGIN{system(\"bsub < /scratch/kannz6/temp/vcfFiles/" + fileName + "\")}' & \n");
	# n = n+1

batchJobsScript.close()
