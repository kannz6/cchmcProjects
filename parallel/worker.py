from time import gmtime, strftime, localtime
import database
import config
import random
import pickle
import shutil
import sys
import os
import re
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

####
#7-13-17
path_to_validated_trios_txt_file = "/scratch/kannz6/loni/vcf_pipeline/done.txt"

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

#################################################################################

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
    
    def _bsub_script_wrap_loni_job(self, W, M, n, J,span, command):
        return "#!/bin/bash\n#BSUB -W %s\n#BSUB -M %s\n#BSUB -n %s\n#BSUB -J %s\n#BSUB -R \"span[%s]\"\n\n%s\n\nexit"%(W, M, n, J,span, command)

    ####################################
    #8-3-17
    def _case_wrap(self,**kwargs):
        if kwargs:
            ret = ""
            if "id" not in kwargs.keys() or "files" not in kwargs.keys() or "finished-file" not in kwargs.keys() or "commands" not in kwargs.keys():
                return "***must pass 'id', 'commands' and 'finished-file' in kwargs dictionary***"
            
            else:
                ##########################################################################
                _test_commands = ["test -e","! test -e","test -s","! test -s","test -n","! test -n"]
                _find_commands = {}; fileSizeList = False;
                # if "search-dir" in kwargs['files']['sizes'].keys():
                #     _searchDir = kwargs['files']['sizes']['search-dir']
                # else:
                _searchDir = "/scratch/kannz6"
                if "sizes" not in kwargs['files'].keys():
                    _find_commands.update({"root":"find {0} -name".format(_searchDir)})
                    _find_commands.update({"output":"find {0} -name".format(self.output_dir)})
                    _find_commands.update({"aligned-sorted":"find {0} -size +2G -name".format(self.output_dir)})
                    _find_commands.update({"aligned-sorted-index":"find {0} -size +2M -name".format(self.output_dir)})
                    _find_commands.update({"sam":"find {0} -size +15G -name".format(_searchDir)})
                elif isinstance(kwargs['files']['sizes'],dict) :
                    fileSizeList = True
                    # _find_command += "find {0} -size {1} -name".format(_searchDir,kwargs['files']['sizes'])

                    [ _find_commands.update({k:"find {0} -size {1} -name".format(_searchDir,v)}) if k not in kwargs['files']['sizes']['search-dirs'].keys() else _find_commands.update({k:"find {0} -size {1} -name".format(kwargs['files']['sizes']['search-dirs'][k],v)}) for k,v in kwargs['files']['sizes'].iteritems() ]

                output_path_sam = "{0}".format(kwargs['files']['sam'])
                output_path_bam = "{0}".format(kwargs['files']['aligned'])
                output_path_sort = "{0}".format(kwargs['files']['aligned-sorted'])
                output_path_sort_index = "{0}".format(kwargs['files']['aligned-sorted-index'])
                output_file_finshed = "{0}/{1}-done.txt".format(self.output_dir,kwargs['id'])
                debug_step_log = "{0}/stepDebug-{1}.txt".format(self.output_dir,self.childs_blinded_id)
                multiKeyBlindedIdFileRename = "{0}/handleMultiKeyBams.sh".format(self.output_dir)
                ##########################################################################
                if len(kwargs['commands']) == 1:
                    ret += "\nc=0\n"
                else:
                    if fileSizeList is True:
                        ret += "\nif ( ( {0} \"{1}\" && {0} \"{2}/$(basename {3})\" && {0} \"{2}/$(basename {4})\" ) || ( ( {7} \"$({5} \"{2}/$(basename {3})\")\" ) && ( {7} \"$({6} \"{2}/$(basename {4})\")\" ) ) )\nthen\n\tc=4\n".format(_test_commands[0],multiKeyBlindedIdFileRename,self.output_dir,output_path_sort,output_path_sort_index,_find_commands['aligned-sorted'],_find_commands['aligned-sorted-index'],_test_commands[0])
                        ret += "\nelif {0} \"$({1} \"$(basename {2})\" -print | egrep -m 1 . | egrep -v $PWD)\"\nthen\n\tpossibleUseableFileDirectory=$(dirname \"$({1} \"$(basename {2} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD)\n\t".format(_test_commands[2],_find_commands['sam'],output_path_sam)

                    else:
                        ret += "if ( ( {0} \"{1}\" && {0} \"{2}/$(basename {3})\" && {0} \"{2}/$(basename {4})\" ) || ( ( {7} \"$({5} \"{2}/$(basename {3})\")\" ) && ( {7} \"$({6} \"{2}/$(basename {4})\" )\" ) ) )\nthen\n\tc=4\n".format(_test_commands[0],multiKeyBlindedIdFileRename,self.output_dir,output_path_sort,output_path_sort_index,_find_commands['aligned-sorted'],_find_commands['aligned-sorted-index'],_test_commands[0])
                        ret += "\nelif {0} \"$({1} \"$(basename {2})\" -print | egrep -m 1 . | egrep -v $PWD)\"\nthen\n\tpossibleUseableFileDirectory=$(dirname \"$({1} \"$(basename {2} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD)\n\t".format(_test_commands[2],_find_commands['sam'],output_path_sam)
                    ret += "if {0} \"$possibleUseableFileDirectory/$(basename {3})\"\n\tthen\n\t\tif {0} \"$possibleUseableFileDirectory/$(basename {4})\"\n\t\tthen\n\t\t\tif ( {5} \"{4}\" && {5} \"{3}\" && {0} \"{6}\" )\n\t\t\tthen\n\t\t\t\tc=3\n\t\t\telif ( {5} \"{3}\" && {5} \"{2}\" && {0} \"{4}\" )\n\t\t\tthen\n\t\t\t\tc=2\n\t\t\telif ( {5} \"{2}\" && {0} \"{3}\" )\n\t\t\tthen\n\t\t\t\tc=1\n\t\t\telse\n\t\t\t\tc=0\n\t\t\tfi\n\t\tfi\n\telse\n\t\tif ( {5} \"{2}\" && {0} \"{3}\" )\n\t\tthen\n\t\t\tc=1\n\t\telse\n\t\t\tc=0\n\t\tfi\n\tfi\nelse\n\tc=0\nfi\n\n".format(_test_commands[2],_find_commands['sam'],output_path_sam,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            output_path_bam,output_path_sort,_test_commands[0],
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            output_path_sort_index)

                ret += "sleep 1\n{0} {1} && echo \"[{3}]- {2} starts at step $c\" > {1} || echo \"[{3}]- {2} starts at step $c\" >> {1}\n\ncase $c in\n\t".format(_test_commands[1],debug_step_log,kwargs['id'],strftime("%H:%M:%S %Y-%m-%d", localtime()))#stopped here 8-11-17

                for n,k in enumerate(kwargs['commands'].keys()):
                    ret += "{0})\n\t".format(n)
                    if n < len(kwargs['commands'].keys()) - 1:
                        ret += "{0}\n\t;;\n\t".format(kwargs['commands'][str(n)])
                    else:
                        ret += "{0}\n\t;;\n".format(kwargs['commands'][str(n)])
            
                ret += "esac\n\nkill %1\n\n"
                
                return ret

        else:
            return "invalid use of _case_wrap method!"         

    def _skip_step_check(self,**kwargs):
        if not kwargs or "check" not in kwargs.keys():
            sys.exit("Bad call to this method! Must at least supply kwarg 'check'...")
        else:
            if "check" in kwargs.keys():
                _test_commands = ["test -e","! test -e","test -s","! test -s","test -n","! test -n"]

                if "file-exists" not in kwargs.keys():
                    _test_command = _test_commands[1]
                else:
                    if kwargs["file-exists"] is True:
                        _test_command = _test_commands[0]
                    else:
                        _test_command = _test_commands[1]

                _find_command = ""; fileSizeList = False;
                if "size" not in kwargs.keys():
                    _find_command += "find /scratch/kannz6 -name"
                else:
                    if "check" in kwargs['size'].keys() and "file-to-copy" in kwargs['size'].keys():
                        fileSizeList = True
                        # _find_command += "find /scratch/kannz6 -size {0} -name".format(kwargs['size'])
                        _find_commands= {}; 
                        [ _find_commands.update({k:"find /scratch/kannz6 -size {0} -name".format(v)}) for k,v in kwargs['size'].iteritems() ]
                    elif "fs" not in kwargs['size'].keys():
                        _find_command += "find /scratch/kannz6 -name"
                    elif "gt" in kwargs['size'].keys() and kwargs['size']['gt'] is True:
                        _find_command += "find /scratch/kannz6 -size +{0} -name".format(kwargs['size']['fs'])
                    elif "gt" in kwargs['size'].keys() and kwargs['size']['gt'] is False:
                        _find_command += "find /scratch/kannz6 -size -{0} -name".format(kwargs['size']['fs'])

                if isinstance(kwargs['check'],list):
                    _ret = "\n\t( "
                    for lc,cfe in enumerate(kwargs['check']):
                        if lc < len(kwargs['check']) - 1:
                            _ret += "( {0} \"{1}\" ) && ".format(_test_command,cfe)
                        else:
                            _ret += "( {0} \"{1}\" )".format(_test_command,cfe)
                    _ret += " ) && ( ( while ( "

                    if "skip-root-directory-check" not in kwargs.keys():
                        kwargs['check'].reverse()
                        for lc,cfe in enumerate(kwargs['check']):
                            if lc < len(kwargs['check']) - 1:
                                _ret += "( {0} \"{1}\" ) && ".format(_test_command,cfe)
                            else:
                                _ret += "( {0} \"{1}\" ) ); do ".format(_test_command,cfe)

                        for lc,cfe in enumerate(kwargs['check']):
                            if lc < len(kwargs['check']) - 1:
                                if fileSizeList is True:
                                    _ret += "existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD) && ( {2} \"$existingFile\" && {3} \"$existingFile\" ) && cp $(dirname $existingFile)/$(basename {1}) {1} || ".format(_find_commands[cfe],cfe,_test_commands[0],_test_commands[2])
                                else:
                                    _ret += "existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD) && ( {2} \"$existingFile\" && {3} \"$existingFile\" ) && cp $(dirname $existingFile)/$(basename {1}) {1} || ".format(_find_command,cfe,_test_commands[0],_test_commands[2])
                            else:
                                if fileSizeList is True:
                                    _ret += "existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD) && ( {2} \"$existingFile\" && {3} \"$existingFile\" ) && cp $(dirname $existingFile)/$(basename {1}) {1} ".format(_find_commands[cfe],cfe,_test_commands[0],_test_commands[2])
                                else:
                                    _ret += "existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD) && ( {2} \"$existingFile\" && {3} \"$existingFile\" ) && cp $(dirname $existingFile)/$(basename {1}) {1} ".format(_find_command,cfe,_test_commands[0],_test_commands[2])
                        _ret +=  "|| break; done ) ||"
                    else:
                        _ret += "&&"
                else:
                    if "if-exists-check" not in kwargs.keys() or ( "if-exists-check" in kwargs.keys() and kwargs['if-exists-check'] is False ) :
                        if "file-to-copy" in kwargs.keys():
                            _ret = "\n\t( {0} \"{1}\" ) ".format(_test_command,kwargs['file-to-copy'])
                        else:
                            _ret = "\n\t( {0} \"{1}\" ) ".format(_test_command,kwargs['check'])
                    else:
                        if kwargs["if-exists-check"]:
                            if "file-to-copy" in kwargs.keys():
                                _ret = "\n\t( {0} \"{1}\" ) ".format(_test_command,kwargs['file-to-copy'])
                            else:
                                _ret = "\n\t( {0} \"{1}\" ) ".format(_test_command,kwargs['check'])
                    if "skip-root-directory-check" not in kwargs.keys():
                        if "size" not in kwargs.keys():
                                                       
                            if "file-to-copy" in kwargs.keys():
                                if fileSizeList is True:
                                    _ret += "&& ( {0} \"$({1} $(basename {2}) -print | egrep -m 1 . | egrep -v $PWD)\" ".format(_test_commands[0],_find_commands['check'],kwargs['check'])
                                    _ret += "&& ( possibleUseableFileDirectory=$(dirname \"$({0} \"$(basename {1} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD) && {2} \"$possibleUseableFileDirectory/$(basename {3})\" ".format(_find_commands['check'],kwargs['check'],_test_commands[0],kwargs['file-to-copy'])

                                else:
                                    _ret += "&& ( {0} \"$({1} $(basename {2}) -print | egrep -m 1 . | egrep -v $PWD)\" ".format(_test_commands[0],_find_command,kwargs['check'])
                                    _ret += "&& ( possibleUseableFileDirectory=$(dirname \"$({0} \"$(basename {1} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD) && {2} \"$possibleUseableFileDirectory/$(basename {3})\" ".format(_find_command,kwargs['check'],_test_commands[0],kwargs['file-to-copy'])
                                _ret += "&& existingFile=$possibleUseableFileDirectory/$(basename {0}) && cp $existingFile {0} ) ||".format(kwargs['file-to-copy'])
                            else:
                                if fileSizeList is True:
                                    _ret += "&& ( existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD| egrep -v $PWD ) && ( {2} \"$existingFile\" && {3} \"$existingFile\") && ".format(_find_command,kwargs['check'],_test_commands[4],_test_commands[2])
                                else:
                                    _ret += "&& ( existingFile=$({0} \"$(basename {1})\" -print | egrep -m 1 . | egrep -v $PWD| egrep -v $PWD ) && ( {2} \"$existingFile\" && {3} \"$existingFile\") && ".format(_find_command,kwargs['check'],_test_commands[4],_test_commands[2])
                                _ret += "cp $existingFile {0} ||".format(kwargs['check'])

                        else:
                            if "file-to-copy" in kwargs.keys():
                                if fileSizeList is True:
                                    _ret += "&& ( {0} \"$({1} $(basename {2}) -print | egrep -m 1 . | egrep -v $PWD)\" ".format(_test_commands[0],_find_commands['check'],kwargs['check'])
                                    _ret += "&& ( possibleUseableFileDirectory=$(dirname \"$({0} \"$(basename {1} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD) && {2} \"$possibleUseableFileDirectory/$(basename {3})\" ".format(_find_commands['check'],kwargs['check'],_test_commands[0],kwargs['file-to-copy'])
                                else:
                                    _ret += "&& ( {0} \"$({1} $(basename {2}) -print | egrep -m 1 . | egrep -v $PWD)\" ".format(_test_commands[0],_find_command,kwargs['check'])
                                    _ret += "&& ( possibleUseableFileDirectory=$(dirname \"$({0} \"$(basename {1} -print | egrep -m 1 . | egrep -v $PWD)\")\" | tac | egrep -m 1 . | egrep -v $PWD) && {2} \"$possibleUseableFileDirectory/$(basename {3})\" ".format(_find_command,kwargs['check'],_test_commands[0],kwargs['file-to-copy'])
                                _ret += "&& existingFile=$possibleUseableFileDirectory/$(basename {0}) && cp $existingFile {0} ) ||".format(kwargs['file-to-copy'])
                            else:
                                _ret += "cp $existingFile {0} ||".format(kwargs['check'])
                    else:
                        _ret += "&&"
                return _ret

    def _in_parallel(self):
        ####
        #4-27-17
        #Testing getting loni id to kill jobs if loni module stopped prematurely
        ####
        loniJob = False
        _kill_job_command = ""
        _job_group_id = ""
        _awk_command = "awk '{print $1}'"
        try:
            if os.path.exists("loniJobId.txt"):
                with open("{0}".format("loniJobId.txt"), "r") as loniIdReader:
                    loniId = loniIdReader.readlines()
                loniIdReader.close()
                loniId = [line.strip() for line in loniId]
                loniJob = True
                _job_group_id = str(loniId[0][3:])
                _initial_parent_status_check = "if ( test -z \"$( bjobs | grep \"{0}\"| {1} )\" ); then bkill $( bjobs | grep \"{2}\"| {1} ); fi".format(loniId[0], _awk_command, _job_group_id)
                _kill_job_command = "while ( ( test -n \"$( bjobs | grep \"{0}\"| {1} )\" ) ); do sleep 1; test -z \"$( bjobs | grep \"{0}\"| {1} )\" && bkill $( bjobs | grep \"{2}\"| {1} ); done &".format(loniId[0], _awk_command, _job_group_id)
                print "[Loni ID: {0}]\ngroupid: {1}".format(str(loniId[0]), str(loniId[0][3:]))

            else:
                print "loniJobId.txt doesn't exist in current directory! Is this a loni job?"
        except Exception as readLoniIdFailed:
            print "[_in_parallel()] Exception thrown when trying to read loniJobId.txt read\nError: {0}\n".format(readLoniIdFailed)
        ####
        bsubCommand = []
        ####################################
        #7-5-17
        # added _used_blinded_ids variable for multi-key-blinded-id
        _used_blinded_ids = []
        for i, pair in enumerate(self.path_pairs):
            commands = "cd {0}\n\t".format(self.home_directory)

            if loniJob:
                commands += "{0}\n\t{1}\n\t".format(_initial_parent_status_check,_kill_job_command)

            blinded_id = pair[0]['id']
            ####################################
            P1_path = os.path.join(config.BASE_PATH, pair[0]['path'][1:] if pair[0]['path'].startswith("/") else pair[0]['path'])
            P2_path = os.path.join(config.BASE_PATH, pair[1]['path'][1:] if pair[0]['path'].startswith("/") else pair[1]['path'])
            ####################################
            #7-5-17
            # added _used_blinded_ids variable for multi-key-blinded-id
            commands5 = commands;
            if len(self.path_pairs) > 3:
                blinded_id += "-{0}".format(i)
                _used_blinded_ids.append(blinded_id)
                chopped_blinded_id = ""
                if len(blinded_id) == 12:
                    chopped_blinded_id = blinded_id[0:10]
                else:
                    chopped_blinded_id = blinded_id[0:7]
                output_path_sam = "{0}/{1}/aligned.{2}.sam".format(self.output_dir,blinded_id,chopped_blinded_id)
                output_path_bam = "{0}/{1}/aligned.{2}.bam".format(self.output_dir,blinded_id,chopped_blinded_id)
                output_path_sort = "{0}/{1}/aligned-sorted.{2}.bam".format(self.output_dir,blinded_id,chopped_blinded_id)            
                output_path_sort_index = "{0}/{1}/aligned-sorted.{2}.bam.bai".format(self.output_dir,blinded_id,chopped_blinded_id)
                commands += "! test -d {0}/{1} && ( mkdir {0}/{1} ) || ( test -n \"$(ls {0}/{1} | grep \"tmp\" | {2})\" && ( rm {0}/{1}/*tmp* && rm {0}/{1}/*sorted* ) || ( ! test -s \"$(ls {0}/{1} | grep \"*bai\" | {2})\" && test -e \"$(ls {0}/{1} | grep \"*sorted*\" | {2})\" && rm \"$(ls {0}/{1} | grep \"*sorted*\" | {2})\") )\n\t".format(self.output_dir,blinded_id,_awk_command)

            #Output Files
            else:
                output_path_sam = "{0}/aligned.{1}.sam".format(self.output_dir, blinded_id)
                output_path_bam = "{0}/aligned.{1}.bam".format(self.output_dir,blinded_id)
                output_path_sort = "{0}/aligned-sorted.{1}.bam".format(self.output_dir,blinded_id)            
                output_path_sort_index = "{0}/aligned-sorted.{1}.bam.bai".format(self.output_dir,blinded_id)

            output_file_finshed = "{0}/{1}-done.txt".format(self.output_dir,blinded_id)
            input_file_vcf = "{0}/{1}.vcf.gz".format(self.output_dir,self.childs_blinded_id)
            completed_mpileup = "{0}/completed-{1}-mpileup.txt".format(self.output_dir,self.childs_blinded_id)
            output_file_filtered_vcf = "{0}/{1}-filtered.vcf".format(self.output_dir,self.childs_blinded_id)
            completed_filtered_vcf = "{0}/completed-{1}-vcf-filter.txt".format(self.output_dir,self.childs_blinded_id)

            ####################################

            ####
            # 3-30-17
            # add dynamic setting of module
            commands2 = commands;

            #add file exist check for file from next step, if it exists, skip this step
            _check_file_exists = self._skip_step_check(**{"check":output_path_bam,"file-to-copy":output_path_sam,"size":{"file-to-copy":"+10G","check":"+2G"}})

            commands += "module load {0}\n\t{1} ( bwa mem -t 8 {2} {3} {4} > {5} ) )\n\t".format(_bwa,_check_file_exists,config.REF_FILE,P1_path,P2_path,output_path_sam)
            commands += "module load {0}\n\t".format(_samtools)

            commands2 += "module load {0}\n\t".format(_samtools)
            commands3 = commands2; commands4 = commands2;

            #add file exist check for file from next step, if it exists, skip this step
            _check_file_exists = self._skip_step_check(**{"check":output_path_sort,"file-to-copy":output_path_bam,"size":{"file-to-copy":"+2G","check":"+2G"}})

            commands += "{0} ( samtools view -Sb {1} > {2} ) )\n\t".format(_check_file_exists,output_path_sam,output_path_bam)
            commands2 += "{0} ( samtools view -Sb {1} > {2} ) )\n\t".format(_check_file_exists,output_path_sam,output_path_bam)

            #comment out 8-10-17
            _check_file_exists = self._skip_step_check(**{"check":output_path_sort_index,"file-to-copy":output_path_sort,"size":{"file-to-copy":"+2G","check":"+2M"}})

            commands += "{0} ( samtools sort -o {1} {2} ) )\n\t".format(_check_file_exists,output_path_sort,output_path_bam)
            commands2 += "{0} ( samtools sort -o {1} {2} ) )\n\t".format(_check_file_exists,output_path_sort,output_path_bam)
            commands3 += "{0} ( samtools sort -o {1} {2} ) )\n\t".format(_check_file_exists,output_path_sort,output_path_bam)

            ###################
            #7-18-17
            #add file exist check for file from next step, if it exists, skip this step

            _check_file_exists = self._skip_step_check(**{"check":output_path_sort_index,"skip-root-directory-check":True})
            commands += "{0} ( samtools index {1} )\n\t".format(_check_file_exists,output_path_sort)
            commands2 += "{0} ( samtools index {1} )\n\t".format(_check_file_exists,output_path_sort)
            commands3 += "{0} ( samtools index {1} )\n\t".format(_check_file_exists,output_path_sort)
            commands4 += "{0} ( samtools index {1} )\n\t".format(_check_file_exists,output_path_sort)#need to remove this extra ')''

            commands += "echo \"complete\" > {0}".format(output_file_finshed)
            commands2 += "echo \"complete\" > {0}".format(output_file_finshed)
            commands3 += "echo \"complete\" > {0}".format(output_file_finshed)
            commands4 += "echo \"complete\" > {0}".format(output_file_finshed)
            commands5 += "echo \"complete\" > {0}".format(output_file_finshed)
            ####3-20-17
            if ( i < len(self.path_pairs) - 1):
                wrap_case_dict = {"id":blinded_id,"finished-file":output_file_finshed,"commands":{"0":commands,"1":commands2,"2":commands3,"3":commands4,"4":commands5}}
                wrap_case_dict.update({"files":{"sam":output_path_sam,"sizes":{"sam":"+15G","aligned-sorted":"+2G","aligned-sorted-index":"+2M","search-dirs":{"aligned-sorted":self.output_dir,"aligned-sorted-index":self.output_dir}},"aligned":output_path_bam,"aligned-sorted":output_path_sort,"aligned-sorted-index":output_path_sort_index}})
                commands = self._case_wrap(**wrap_case_dict)

            elif( i == len(self.path_pairs) - 1):
                child_id = self.childs_blinded_id
                mom_id = self.childs_blinded_id + '-01'
                dad_id = self.childs_blinded_id + '-02'
                ####################################
                #7-5-17
                # added _used_blinded_ids variable and if-elif logic for handling multi-key-blinded-id 
                if len(self.path_pairs) == 3:
                    commands += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,child_id, mom_id, dad_id)
                    commands2 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,child_id, mom_id, dad_id)
                    commands3 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,child_id, mom_id, dad_id)
                    commands4 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,child_id, mom_id, dad_id)
                    commands5 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" && test -e \"{0}/{2}-done.txt\" && test -e \"{0}/{3}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,child_id, mom_id, dad_id)

                elif len(self.path_pairs) > 3:

                    for j,b_id_multi_key in enumerate(_used_blinded_ids):
                        if j == 0:
                            commands += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands2 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands3 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands4 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands5 += "\n\twhile ! (test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                        elif j < (len(_used_blinded_ids) - 1) :
                            commands += "&& test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands2 += "&& test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands3 += "&& test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands4 += "&& test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                            commands5 += "&& test -e \"{0}/{1}-done.txt\" ".format(self.output_dir,b_id_multi_key)
                        else:
                            commands += "&& test -e \"{0}/{1}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,b_id_multi_key)
                            commands2 += "&& test -e \"{0}/{1}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,b_id_multi_key)
                            commands3 += "&& test -e \"{0}/{1}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,b_id_multi_key)
                            commands4 += "&& test -e \"{0}/{1}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,b_id_multi_key)
                            commands5 += "&& test -e \"{0}/{1}-done.txt\"); do sleep 180; done;\n\t".format(self.output_dir,b_id_multi_key)
                            _rename_multi_key_python_script = "multiKeyBlindedIdFileRename.py"
                            _comleted_multi_key_handle_file = "{0}/completed-handleMultiKeyBams.txt".format(self.output_dir)
                            commands += "cp {0} {1}\n\t".format(_rename_multi_key_python_script,self.output_dir)
                            commands += "cd {0}\n\t".format(self.output_dir)
                            commands += "./{0}\n\t".format(_rename_multi_key_python_script)
                            commands += "bsub < handleMultiKeyBams.sh\n\t"
                            commands += "cd ..\n\t\n\twhile ! (test -e \"{0}\"); do sleep 36; done;\n\n\t".format(_comleted_multi_key_handle_file)
                            commands2 += "cp {0} {1}\n\t".format(_rename_multi_key_python_script,self.output_dir)
                            commands2 += "cd {0}\n\t".format(self.output_dir)
                            commands2 += "./{0}\n\t".format(_rename_multi_key_python_script)
                            commands2 += "bsub < handleMultiKeyBams.sh\n\t"
                            commands2 += "cd ..\n\t\n\twhile ! (test -e \"{0}\"); do sleep 36; done;\n\n\t".format(_comleted_multi_key_handle_file)
                            commands3 += "cp {0} {1}\n\t".format(_rename_multi_key_python_script,self.output_dir)
                            commands3 += "cd {0}\n\t".format(self.output_dir)
                            commands3 += "./{0}\n\t".format(_rename_multi_key_python_script)
                            commands3 += "bsub < handleMultiKeyBams.sh\n\t"
                            commands3 += "cd ..\n\t\n\twhile ! (test -e \"{0}\"); do sleep 36; done;\n\n\t".format(_comleted_multi_key_handle_file)
                            commands4 += "cp {0} {1}\n\t".format(_rename_multi_key_python_script,self.output_dir)
                            commands4 += "cd {0}\n\t".format(self.output_dir)
                            commands4 += "./{0}\n\t".format(_rename_multi_key_python_script)
                            commands4 += "bsub < handleMultiKeyBams.sh\n\t"
                            commands4 += "cd ..\n\t\n\twhile ! (test -e \"{0}\"); do sleep 36; done;\n\n\t".format(_comleted_multi_key_handle_file)
                            # commands5 += "cp {0} {1}\n\t".format(_rename_multi_key_python_script,self.output_dir)
                            # commands5 += "cd {0}\n\t".format(self.output_dir)
                            # commands5 += "./{0}\n\t".format(_rename_multi_key_python_script)
                            # commands5 += "bsub < handleMultiKeyBams.sh\n\t"
                            # commands5 += "cd ..\n\t\n\twhile ! (test -e \"{0}\"); do sleep 36; done;\n\n\t".format(_comleted_multi_key_handle_file)
                ############################################################################
                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "module load {0}\n\tmodule load {1}\n\t".format(_samtools, _bcftools)
                commands2 += "module load {0}\n\tmodule load {1}\n\t".format(_samtools, _bcftools)
                commands3 += "module load {0}\n\tmodule load {1}\n\t".format(_samtools, _bcftools)
                commands4 += "module load {0}\n\tmodule load {1}\n\t".format(_samtools, _bcftools)
                commands5 += "module load {0}\n\tmodule load {1}\n\t".format(_samtools, _bcftools)
                ############################################################################
   
                _check_file_exists = self._skip_step_check(**{"check":completed_filtered_vcf,"file-to-copy":output_file_filtered_vcf,"size":{"file-to-copy":"+19M","check":"+36"}})

                commands += "{0} ( samtools mpileup -t AD -C50 -uf {1} {2}/aligned-sorted.{3}.bam {2}/aligned-sorted.{4}.bam {2}/aligned-sorted.{5}.bam".format(_check_file_exists,config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands2 += "{0} ( samtools mpileup -t AD -C50 -uf {1} {2}/aligned-sorted.{3}.bam {2}/aligned-sorted.{4}.bam {2}/aligned-sorted.{5}.bam".format(_check_file_exists,config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands3 += "{0} ( samtools mpileup -t AD -C50 -uf {1} {2}/aligned-sorted.{3}.bam {2}/aligned-sorted.{4}.bam {2}/aligned-sorted.{5}.bam".format(_check_file_exists,config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands4 += "{0} ( samtools mpileup -t AD -C50 -uf {1} {2}/aligned-sorted.{3}.bam {2}/aligned-sorted.{4}.bam {2}/aligned-sorted.{5}.bam".format(_check_file_exists,config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands5 += "{0} ( samtools mpileup -t AD -C50 -uf {1} {2}/aligned-sorted.{3}.bam {2}/aligned-sorted.{4}.bam {2}/aligned-sorted.{5}.bam".format(_check_file_exists,config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)

                commands += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} ) ) || ".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands2 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} ) ) || ".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands3 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} ) ) || ".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands4 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} ) ) || ".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands5 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} ) ) || ".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)

                check_filtered_file_size_command = "$(find {0} -size +15M -name \"{1}-filtered.vcf\" -print)".format(self.output_dir,self.childs_blinded_id)
                commands += "( ( ! test -s \"{0}\" ) && samtools mpileup -t AD -C50 -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands2 += "( ( ! test -s \"{0}\" ) && samtools mpileup -t AD -C50 -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands3 += "( ( ! test -s \"{0}\" ) && samtools mpileup -t AD -C50 -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands4 += "( ( ! test -s \"{0}\" ) && samtools mpileup -t AD -C50 -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)
                commands5 += "( ( ! test -s \"{0}\" ) && samtools mpileup -t AD -C50 -uf {0} {1}/aligned-sorted.{2}.bam {1}/aligned-sorted.{3}.bam {1}/aligned-sorted.{4}.bam".format(config.REF_FILE,self.output_dir,child_id,mom_id,dad_id)

                commands += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} || echo \"{1}'s vcf file has been filtered.\" > {2} )\n\t".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands2 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} || echo \"{1}'s vcf file has been filtered.\" > {2} )\n\t".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands3 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} || echo \"{1}'s vcf file has been filtered.\" > {2} )\n\t".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands4 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} || echo \"{1}'s vcf file has been filtered.\" > {2} )\n\t".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)
                commands5 += " | bcftools call -mv -Ou | bcftools filter -g3 -G10 -e'AN<5 || %QUAL<10 || %MAX(AD[1])<=3 || %MAX(AD[1])/%MAX(DP)<=0.144' | bcftools view -m2 -M2 -v snps > {0} && echo \"{1}'s vcf file has been filtered.\" > {2} || echo \"{1}'s vcf file has been filtered.\" > {2} )\n\t".format(output_file_filtered_vcf,self.childs_blinded_id,completed_filtered_vcf)

                local_error_output = "{0}/{1}-error.txt".format(self.output_dir,self.childs_blinded_id)
                

                commands += "test -s \"{0}\" && rm $(find -type f -name \"*sam\") || ( echo \"Issue during mpileup step. Filtered vcf file is corrupted or doesn't exist!\" > {1} && kill %1 )\n\t".format(check_filtered_file_size_command,local_error_output)
                commands2 += "test -s \"{0}\" && rm $(find -type f -name \"*sam\") || ( echo \"Issue during mpileup step. Filtered vcf file is corrupted or doesn't exist\" > {1} && kill %1 )\n\t".format(check_filtered_file_size_command,local_error_output)
                commands3 += "test -s \"{0}\" && rm $(find -type f -name \"*sam\") || ( echo \"Issue during mpileup step. Filtered vcf file is corrupted or doesn't exist\" > {1} && kill %1 )\n\t".format(check_filtered_file_size_command,local_error_output)
                commands4 += "test -s \"{0}\" && rm $(find -type f -name \"*sam\") || ( echo \"Issue during mpileup step. Filtered vcf file is corrupted or doesn't exist\" > {1} && kill %1 )\n\t".format(check_filtered_file_size_command,local_error_output)
                commands5 += "test -s \"{0}\" && rm $(find -type f -name \"*sam\") || ( echo \"Issue during mpileup step. Filtered vcf file is corrupted or doesn't exist\" > {1} && kill %1 )\n\t".format(check_filtered_file_size_command,local_error_output)

                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "module load {0}\n\t".format(_plink)
                commands2 += "module load {0}\n\t".format(_plink)
                commands3 += "module load {0}\n\t".format(_plink)
                commands4 += "module load {0}\n\t".format(_plink)
                commands5 += "module load {0}\n\t".format(_plink)
                ####
                output_file_plink = "{0}/{1}.plink".format(self.output_dir,self.childs_blinded_id)
                commands += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n\t".format(output_file_filtered_vcf,output_file_plink)
                commands2 += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n\t".format(output_file_filtered_vcf,output_file_plink)
                commands3 += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n\t".format(output_file_filtered_vcf,output_file_plink)
                commands4 += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n\t".format(output_file_filtered_vcf,output_file_plink)
                commands5 += "plink --allow-extra-chr --vcf {0} --make-bed --out {1}\n\t".format(output_file_filtered_vcf,output_file_plink)

                ##############
                #4-7-2017
                #Add Andys changes where we update the plink file using the blind-ids given in *plink.fam file
                ########
                fam_file_fixer = "update_fam.py"
                plink_fam_file = "{0}.plink.fam".format(self.childs_blinded_id)
                commands += "cp {0} {1}\n\tcd {1}\n\t./{0} {2}\n\trm {0}\n\tcd ..\n\t".format(fam_file_fixer,self.output_dir,plink_fam_file)
                commands2 += "cp {0} {1}\n\tcd {1}\n\t./{0} {2}\n\trm {0}\n\tcd ..\n\t".format(fam_file_fixer,self.output_dir,plink_fam_file)
                commands3 += "cp {0} {1}\n\tcd {1}\n\t./{0} {2}\n\trm {0}\n\tcd ..\n\t".format(fam_file_fixer,self.output_dir,plink_fam_file)
                commands4 += "cp {0} {1}\n\tcd {1}\n\t./{0} {2}\n\trm {0}\n\tcd ..\n\t".format(fam_file_fixer,self.output_dir,plink_fam_file)
                commands5 += "cp {0} {1}\n\tcd {1}\n\t./{0} {2}\n\trm {0}\n\tcd ..\n\t".format(fam_file_fixer,self.output_dir,plink_fam_file)

                tmp_fam_file = "{0}/tmp_fam.txt".format(self.output_dir)
                commands += "cp {0} {1}/{2}\n\t".format(tmp_fam_file,self.output_dir,plink_fam_file)
                commands2 += "cp {0} {1}/{2}\n\t".format(tmp_fam_file,self.output_dir,plink_fam_file)
                commands3 += "cp {0} {1}/{2}\n\t".format(tmp_fam_file,self.output_dir,plink_fam_file)
                commands4 += "cp {0} {1}/{2}\n\t".format(tmp_fam_file,self.output_dir,plink_fam_file)
                commands5 += "cp {0} {1}/{2}\n\t".format(tmp_fam_file,self.output_dir,plink_fam_file)

                ####
                # 3-30-17
                # add dynamic setting of module
                commands += "module load {0}\n\t".format(_king)
                commands2 += "module load {0}\n\t".format(_king)
                commands3 += "module load {0}\n\t".format(_king)
                commands4 += "module load {0}\n\t".format(_king)
                commands5 += "module load {0}\n\t".format(_king)
                ####
                output_file_king = "{0}/{1}.king".format(self.output_dir,self.childs_blinded_id)
                commands += "king -b {0}.bed --kinship --prefix {1}\n\t".format(output_file_plink,output_file_king)
                commands2 += "king -b {0}.bed --kinship --prefix {1}\n\t".format(output_file_plink,output_file_king)
                commands3 += "king -b {0}.bed --kinship --prefix {1}\n\t".format(output_file_plink,output_file_king)
                commands4 += "king -b {0}.bed --kinship --prefix {1}\n\t".format(output_file_plink,output_file_king)
                commands5 += "king -b {0}.bed --kinship --prefix {1}\n\t".format(output_file_plink,output_file_king)

                commands += "#cleanup directory for space management\n\ttest -s \"{0}\" && rm {1}/*plink*\n\ttest -s \"{0}\" && rm {1}/*.sam\n\trm {1}/*TMP*\n\trm {1}/*kingX*\n\t".format(tmp_fam_file,self.output_dir)         
                commands2 += "#cleanup directory for space management\n\t\n\ttest -s \"{0}\" && rm {1}/*plink*\n\ttest -s \"{0}\" && rm {1}/*.sam\n\trm {1}/*TMP*\n\trm {1}/*kingX*\n\t".format(tmp_fam_file,self.output_dir)         
                commands3 += "#cleanup directory for space management\n\t\n\ttest -s \"{0}\" && rm {1}/*plink*\n\ttest -s \"{0}\" && rm {1}/*.sam\n\trm {1}/*TMP*\n\trm {1}/*kingX*\n\t".format(tmp_fam_file,self.output_dir)         
                commands4 += "#cleanup directory for space management\n\t\n\ttest -s \"{0}\" && rm {1}/*plink*\n\ttest -s \"{0}\" && rm {1}/*.sam\n\trm {1}/*TMP*\n\trm {1}/*kingX*\n\t".format(tmp_fam_file,self.output_dir)
                commands5 += "#cleanup directory for space management\n\t\n\ttest -s \"{0}\" && rm {1}/*plink*\n\ttest -s \"{0}\" && rm {1}/*.sam\n\trm {1}/*TMP*\n\trm {1}/*kingX*\n\t".format(tmp_fam_file,self.output_dir)         

                
                output_job_finshed = "{0}/{1}-trio-validation-complete.txt".format(self.output_dir,self.childs_blinded_id)
                # commands += "echo \"complete\" > {0}\nkill %1".format(output_job_finshed)#7-13-17 comment out for test line 335
                commands += "test -s \"{0}.kin\" && ( echo \"complete\" > {1} && echo \"{2}\" >> {3} )\n\t".format(output_file_king,output_job_finshed,self.childs_blinded_id,path_to_validated_trios_txt_file)
                commands2 += "test -s \"{0}.kin\" && ( echo \"complete\" > {1} && echo \"{2}\" >> {3} )\n\t".format(output_file_king,output_job_finshed,self.childs_blinded_id,path_to_validated_trios_txt_file)
                commands3 += "test -s \"{0}.kin\" && ( echo \"complete\" > {1} && echo \"{2}\" >> {3} )\n\t".format(output_file_king,output_job_finshed,self.childs_blinded_id,path_to_validated_trios_txt_file)
                commands4 += "test -s \"{0}.kin\" && ( echo \"complete\" > {1} && echo \"{2}\" >> {3} )\n\t".format(output_file_king,output_job_finshed,self.childs_blinded_id,path_to_validated_trios_txt_file)
                commands5 += "test -s \"{0}.kin\" && ( echo \"complete\" > {1} && echo \"{2}\" >> {3} )\n\t".format(output_file_king,output_job_finshed,self.childs_blinded_id,path_to_validated_trios_txt_file)
                

                wrap_case_dict = {"id":blinded_id,"finished-file":output_file_finshed,"commands":{"0":commands,"1":commands2,"2":commands3,"3":commands4,"4":commands5}}
                wrap_case_dict.update({"files":{"sam":output_path_sam,"sizes":{"sam":"+15G","aligned-sorted":"+2G","aligned-sorted-index":"+2M","search-dirs":{"aligned-sorted":self.output_dir,"aligned-sorted-index":self.output_dir}},"aligned":output_path_bam,"aligned-sorted":output_path_sort,"aligned-sorted-index":output_path_sort_index}})
                commands = self._case_wrap(**wrap_case_dict)
            ####
            #5-2-17
            #add job name of loni job for killing
            if loniJob:
                script = self._bsub_script_wrap_loni_job(str(_wall),str(_memory),str(_cores), _job_group_id ,("ptile=%s"%(_cores)),commands)
            else:
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
