import sys
import os

import config
import database

sys.path.insert(0, config.LSF_PATH)
import LSF

import pickle

if __name__ == "__main__":

    print "***main:"

    print "Loading fastq.gz file paths from postgres..."
    # load both the fastq.gz paths and the blinded ids
    ids_and_paths = database.get_trio_file_paths()

    print "collecting all ids..."
    # get a list of all blinded ids
    all_ids = map(lambda x: x[0], ids_and_paths)
    # print "all ids: %s"% all_ids

    print "collecting children..."
    # get a set of all blinded ids that are not a father or mother, according to the id convention.
    children = set(filter(lambda x: not x.endswith('-01') and not x.endswith('-02'), all_ids))
    # print "children: %s"%children

    print "collecting moms..."
    # get a set of all blinded ids of mothers, according to the id convention.
    moms = set(filter(lambda x: x.endswith('-01'), all_ids))
    # print "moms: %s"% moms

    print "collecting dads..."
    # get a set of all blinded ids of fathers, according to the id convention.
    dads = set(filter(lambda x: x.endswith('-02'), all_ids))
    # print "dads: %s"% dads

    print "collecting children of trios..."
    # select each child id such that the mom and the dad exist (its a trio).
    trio_children = filter(lambda x: x + '-01' in moms and x + '-02' in dads, children)

    print "collecting trios..."
    # get the ids and file paths for the known trios.
    trios = [filter(lambda x: x[0] in (ID,ID+'-01',ID+'-02'),ids_and_paths) for ID in trio_children]
    # print trios
    # print len(trios)
    # serialize the records so that they can be written to a file.
    trios = [pickle.dumps(trio) for trio in trios]

    print "Running cluster job..."
    lsf_job = LSF.LSFjob(worker_script='worker.py',
                         operands=trios,
                         node_count=len(trios),
                         RAM=5000, # MB
                         run_ceiling_hours=24,
                         run_ceiling_minutes=0,
                         delim=',') # delimit records in the input files using commas because pickle uses newline characters.
    lsf_job.show_params()
    lsf_job.run()
