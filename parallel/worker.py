from time import gmtime, strftime, localtime
import database
import config
import random
import pickle
import shutil
import sys
import os
####3-17-17
from os import getpid, path
import versionSettings

#####
# David Beatrice
# added levenshtein distance and usage
# levenshtein distance between two strings
# (taken from wikipedia)
def lev(s1, s2):
    if len(s1) < len(s2):
        return lev(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]
# 

####
# 3-24-17
number_of_operands = 0
desired_number_of_cores = 10
desired_memory = 18000

#####
# 3-30-17
# Add dynamic version settings
_wall = versionSettings.wall
_memory = versionSettings.memory
_bwa = "bwa/{0}".format(versionSettings.bwa)
_samtools = "samtools/{0}".format(versionSettings.samtools)
_bcftools = "bcftools/{0}".format(versionSettings.bcftools)
_plink = "plink/{0}".format(versionSettings.plink)
_king = "king/{0}".format(versionSettings.king)
####

####
# 3-29-17
# Add setCores
def setCores(v):
    # print "inside setCores\nv: %s"% v
    if ( v < 25):
        return 14
    elif( v >= 25 and v < 50 ):
        return 12
    elif( v >= 50 and v < 75):
        return 10
    elif(v >= 75 ):
        return 8

x = 0
########

class Trio:
    def __init__(self, path_records):
        #set the name of the executing directory 
        self.home_directory = os.getcwd()
        # each record consists of a path to a fastq file, and its respective blinded id
        self.path_records = self._to_dicts(path_records)
        # this list will hold P1/P2 pairs of fastq files
        self.path_pairs = []
        # where to print output files
        self.output_dir = ''
        # save the child's blinded id here - the parent's ids can be inferred from it
        self.childs_blinded_id = ''
        # where to print error messages
        self.error_log_file = ''
        # register pre-processing errors here so that we know not to proceed with the pipeline
        self.errors = False

        self._get_childs_blinded_id()
        self._create_output_directory()
        self._set_log_file()
        # make sure that we have the minimum requred files to run the pipeline successfully
        self._check_minimal_file_requirements()
        if not self.errors:
            # pair up the P1 and P2 files
            self._pair_paths()

    def _log_error(self,message):
        path = os.path.join(self.output_dir,self.error_log_file)
        with open(path,'a') as f:
            f.write("(%s) %s\n"%(strftime("%Y-%m-%d %H:%M:%S", gmtime()),message))

    def _set_log_file(self):
        self.error_log_file = "error_log.%s.txt"%(self.childs_blinded_id)

    def _check_minimal_file_requirements(self):
        # Make sure that we have at least one P1 and one P2 fastq file for each member of the trio.
        # Also, check that all file paths point to real files.

        # get parent's blinded ids
        mom_id = self.childs_blinded_id + '-01'
        dad_id = self.childs_blinded_id + '-02'

        # seperate the P1 files from the P2 files
        P1s = filter(lambda x: x['path'].endswith('P1.fastq.gz'), self.path_records)
        P2s = filter(lambda x: x['path'].endswith('P2.fastq.gz'), self.path_records)

        # get the ids for each group
        P1_ids = set(map(lambda x: x['id'], P1s))
        P2_ids = set(map(lambda x: x['id'], P2s))

        # make sure each person has both a P1 and a P2 fastq file
        for _id in [self.childs_blinded_id,mom_id,dad_id]:
            if not (_id in P1_ids and _id in P2_ids):
                self._log_error("Missing P1 or P2 fastq.gz file for blinded id, %s."%(_id))
                self.errors = True

        # make sure all the files exist
        for record in self.path_records:
            path = os.path.join(config.BASE_PATH, record['path'][1:] if record['path'].startswith("/") else record['path'])
            if not os.path.isfile(path):
                self._log_error("Missing file: %s"%(path))
                self.errors = True

    def _to_dicts(self,path_records):
        # re-shape the input data
        return map(lambda x: {'id':x[0], 'path':x[1]}, path_records)

    def _get_childs_blinded_id(self):
        # obtain the child's blinded id for labeling the output
        _id = self.path_records[0]['id']
        if _id.endswith('-01') or _id.endswith('-02'):
            self.childs_blinded_id = _id[:-3]
        else:
            self.childs_blinded_id = _id

    def _create_output_directory(self):
        # Create an output directory named after the child's blind id.
        self.output_dir = self.childs_blinded_id + "_output"
        # create the directory if it DNE
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

    def _pair_paths(self):
        # Sometimes a given trio member will have multiple P1.fastq and P2.fastq files.
        # This method pairs P1s and P2s up according to how likely they were to have
        # been generated in the same batch - using file name similarity as a proxy.

        # seperate the p1s from the p2s because we're going to match them.
        P1s = filter(lambda x: x['path'].endswith('P1.fastq.gz'), self.path_records)
        P2s = filter(lambda x: x['path'].endswith('P2.fastq.gz'), self.path_records)
        # for each p1, find the best matching p2
        for p1_index, p1 in enumerate(P1s):
            # Only consider match candidates with the same blinded ids
            # (meaning, only match fathers to fathers, mothers to mothers, etc...).
            P2_id_matches = filter(lambda x: x['id'] == p1['id'],P2s)
            # Find a p2 to pair with the p1 such that their file names are the best match.
            distances = map(lambda x: lev(x,p1), P2_id_matches)
            # get the index of the matching p2
            if len(distances) == 0:
                print P2s
            p2_index = distances.index(min(distances))
            # add the pair to our list.
            self.path_pairs.append([P1s[p1_index],P2_id_matches[p2_index]])

    def _run(self,commands):
        os.system(" && ".join(commands))
        # print out the command for debugging purposes
        print "\n".join(commands)

    def _bsub_script_wrap(self, W, M, n, span, command):
        return "#!/bin/bash\n#BSUB -W %s\n#BSUB -M %s\n#BSUB -n %s\n#BSUB -R \"span[%s]\"\n\n%s\n\nexit"%(W, M, n, span, command)

    def _in_parallel(self):
        ####
        #4-27-17
        #Testing getting loni id to kill jobs if loni module stopped prematurely
        ####
        # loniJob = False
        # _kill_job_command = ""
        # try:
        #     if os.path.exists("loniJobId.txt"):
        #         with open("{0}".format("loniJobId.txt"), "r") as loniIdReader:
        #             loniId = loniIdReader.readlines()
        #         loniIdReader.close()
        #         loniId = [line.strip() for line in loniId]
        #         loniJob = True
        #         _awk_command = "awk '{print $1}'"
        #         _kill_job_command = "while ( ( test -n \"$( bjobs | grep \"{0}\"| {1} )\" ) ); do sleep 1; test -z \"$( bjobs | grep \"{0}\"| {1} )\" && exit; done &".format(loniId[0], _awk_command)
        #         print "[Loni ID: {0}]\n".format(str(loniId[0]))
        #     else:
        #         print "loniJobId.txt doesn't exist in current directory! Is this a loni job?"
        # except Exception as readLoniIdFailed:
        #     print "[_in_parallel()] Exception thrown when trying to read loniJobId.txt read\nError: {0}\n".format(readLoniIdFailed)
        # ####
        bsubCommand = []
        for i, pair in enumerate(self.path_pairs):
            commands = ""

            blinded_id = pair[0]['id']
            P1_path = os.path.join(config.BASE_PATH, pair[0]['path'][1:] if pair[0]['path'].startswith("/") else pair[0]['path'])
            P2_path = os.path.join(config.BASE_PATH, pair[1]['path'][1:] if pair[0]['path'].startswith("/") else pair[1]['path'])
            output_path_sam = "{0}/aligned.{1}.sam".format(self.output_dir, blinded_id)
            commands += "cd {0}\n".format(self.home_directory)        
            ####
            # 3-30-17
            # add dynamic setting of module
            # if loniJob:
            #     _backroundJobId = "$!"
            #     _backroundJobFile = "backroundID.txt"
            #     commands += "{0}\n".format(_kill_job_command)
            #     if i == 0:
            #         commands += "echo \"{0}\" > {1}\n\n".format( _backroundJobId,_backroundJobFile)
            #     else:
            #         commands += "echo \"{0}\" >> {1}\n\n".format( _backroundJobId,_backroundJobFile)
            commands += "module load {0}\nbwa mem -t 8 {1} {2} {3} > {4}\n".format(_bwa,config.REF_FILE,P1_path,P2_path,output_path_sam)
            commands += "module load {0}\n".format(_samtools)
            ####
            output_path_bam = "{0}/aligned.{1}.bam".format(self.output_dir,blinded_id)
            commands += "samtools view -Sb {0} > {1}\n".format(output_path_sam,output_path_bam)

            ####3-20-17
            commands += "rm {0}\n".format(output_path_sam)
            ####
            output_path_sort = "{0}/aligned-sorted.{1}.bam".format(self.output_dir,blinded_id)
            commands += "samtools sort -o {0} {1}\n".format(output_path_sort,output_path_bam)

            commands += "samtools index {0}\n".format(output_path_sort)

            output_file_finshed = "{0}/{1}-done.txt".format(self.output_dir,blinded_id)
            commands += "echo \"complete\" > {0}".format(output_file_finshed)
            
            ####3-20-17
            if( i == len(self.path_pairs)-1):
                child_id = self.childs_blinded_id
                mom_id = self.childs_blinded_id + '-01'
                dad_id = self.childs_blinded_id + '-02'
                commands += "\nwhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n".format(self.output_dir,child_id, mom_id, dad_id)
                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "#module load {0}\nmodule load {1}\n".format(_samtools, _bcftools)
                ####
                commands += "samtools mpileup -t AD -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands += " | bcftools call -mv -Oz > {0}/{1}.vcf.gz\n".format(self.output_dir,self.childs_blinded_id)

                input_file_vcf = "{0}/{1}.vcf.gz".format(self.output_dir,self.childs_blinded_id)
                output_file_filtered_vcf = "{0}/{1}-filtered.vcf".format(self.output_dir,self.childs_blinded_id)
                # 4-12-17
                #added pipe to bcftools view command
                commands += "bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' {0} | bcftools view -m2 -M2 -v snps > {1}\n".format(input_file_vcf, output_file_filtered_vcf)
                # commands += "bcftools view -m2 -M2 -v snps {0} > {1}\n".format(tmp_filtered_vcf, output_file_filtered_vcf)
                commands += "rm {0}/{1}.vcf.gz\n".format(self.output_dir,self.childs_blinded_id)

                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "module load {0}\n".format(_plink)
                ####
                output_file_plink = "{0}/{1}.plink".format(self.output_dir,self.childs_blinded_id)
                commands += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n".format(output_file_filtered_vcf,output_file_plink)

                ##############
                #4-7-2017
                #Add Andys changes where we update the plink file using the blind-ids given in *plink.fam file
                ########
                fam_file_fixer = "update_fam.py"
                plink_fam_file = "{0}.plink.fam".format(self.childs_blinded_id)
                commands += "cp {0} {1}\ncd {1}\n./{0} {2}\nrm {0}\ncd ..\n".format(fam_file_fixer,self.output_dir,plink_fam_file)
                tmp_fam_file = "{0}/tmp_fam.txt".format(self.output_dir)
                commands += "cp {0} {1}/{2}\n".format(tmp_fam_file,self.output_dir,plink_fam_file)
                # comment out echoing of .fam file fixing 4-7-17
                # commands += "echo \"{1} {1} {1}-01 {1}-02 0 -9\n{1} {1}-01 0 0 2 -9\n{1} {1}-02 0 0 1 -9\" > {0}/{1}.plink.fam\n".format(self.output_dir,self.childs_blinded_id)

                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "module load {0}\n".format(_king)
                ####
                output_file_king = "{0}/{1}.king".format(self.output_dir,self.childs_blinded_id)
                commands += "king -b {0}.bed --kinship --prefix {1}\n".format(output_file_plink,output_file_king)
                commands += "#cleanup directory for space management\nrm {0}/*plink*\nrm {0}/*bam\nrm {0}/*.bai\nrm {0}/*.vcf.gz\nrm {0}/*TMP*\nrm {0}/*kingX*\n".format(self.output_dir)         
                
                output_job_finshed = "{0}/{1}-trio-validation-complete.txt".format(self.output_dir,blinded_id)
                commands += "echo \"complete\" > {0}\n\n".format(output_job_finshed)
            ####
            # 3-29-17
            script = self._bsub_script_wrap(str(_wall),str(_memory),str(_cores),("ptile=%s"%(_cores)),commands)
            ####
            scriptName = "{0}/parallelize_{1}.sh".format(self.output_dir,blinded_id)

            with open( scriptName, 'w' ) as batchScript:
                batchScript.write(script)
            batchScript.close()
            
            bsubCommand.append("bsub < {0}/{1}".format(self.home_directory,scriptName))

        self._run(bsubCommand)

    def verify(self):
        if not self.errors:
            self._in_parallel()

def parse_command_line_arg():
    if len(sys.argv) < 2:
        raise Exception("This script requires an input file path argument.")
    input_file = sys.argv[1]

    #####
    #3-24-17
    global number_of_operands; global x; global _cores;
    number_of_operands = sys.argv[2]
    # 3-29-17
    # add set cores call depending on number of operands
    x = ((int(number_of_operands) + 1) * desired_number_of_cores)
    _cores = setCores(x)
    ######

    # check if its a valid file
    if not os.path.isfile(input_file):
        raise Exception("%s is not a valid file" % (input_file))
    return input_file

if __name__ == "__main__":
    # get the input file
    input_file = parse_command_line_arg()
    # get the paths and blinded ids
    with open(input_file,'r') as f:
        fastq_paths = [pickle.loads(line) for line in f.read().split(',')]

    
    for paths in fastq_paths:
        trio = Trio(paths)
        #####
        #3-24-17
        print "# of operands: {0}".format(number_of_operands)
        print "ptile= {0}".format(_cores)
        print "cores= {0}".format(_cores)
        print "[module versions]\nbwa: {0}\nsamtools: {1}\nbcftools: {2}\nplink: {3}\nking: {4}\n".format(_bwa,_samtools,_bcftools,_plink,_king)
        #####
        #3-24-17
        trio.verify()