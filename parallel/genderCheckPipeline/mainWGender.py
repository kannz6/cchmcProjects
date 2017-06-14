import config
import database
import math
import os
import pickle
import sys

sys.path.insert(0, config.LSF_PATH)
import LSF

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
    triosUsedChecker = open("triosStructure.txt", "w+")
    triosUsedChecker.write(str(trios))

    # print trios
    # print len(trios)
    # serialize the records so that they can be written to a file.
    # trios = [pickle.dumps(trio) for trio in trios]#comment out 5-25-17
    loni_job = False
    try:
        with open("{0}".format("loniJobId.txt"), "r") as loniIdReader:
            loniId = loniIdReader.readlines()
        loniIdReader.close()
        loniId = [line.strip() for line in loniId]
        _job_group_id = str(loniId[0][3:])
        loni_job = True
    except Exception as _no_loni_file_found:
        pass

    _max_node_count = 70
    print "Running cluster job..."
    if len(trios) > _max_node_count:
            _node_count = 3
            triosUsedChecker.write("\n\ninside if len(trios) > 70".format(len(trios)))
            _loops = math.ceil(len(trios)/_node_count + 1); x = 0; y = _node_count
            _count = 0
            n = 0
            triosUsedChecker.write("\n\n[lsf_jobs]\n")
            while _count != _loops:
                if _count == 0:
                    _trios_slice = trios[x:y]
                else:
                    try:
                        if y+1 <= (len(trios) - 1):
                            _trios_slice = trios[x:y+1]
                        else:
                            _trios_slice = trios[x:y]
                    except Exception as _bad_trios_index:
                        break
                try:
                    # _ids = [ v[0][0] for v in _trios_slice if "-02" not in v and "-01" not in v ]#works
                    _ids = [ [ _id[0] for c,_id in enumerate(v) if c % 2 == 0] for v in _trios_slice ]
                    _trio = _trios_slice[0][0]
                except Exception as _bad_index:
                    break
                triosUsedChecker.write("\n---------------------\ncount:{0}\n\nslice: {1}\n\n\n_trio: {2}\n\n\n_ids: {3}\n\n{4}".format(_count,str(_trios_slice),_trio,_ids,_ids[len(_ids)-1][0]))

                depp = "{0}_output/{0}-trio-validation-complete.txt".format(_trio[0])

                _trios_slice = [pickle.dumps(trio) for trio in _trios_slice]
                lsf_job = LSF.LSFjob(worker_script='workerWGender.py',
                                 operands=_trios_slice,
                                 node_count=len(_trios_slice),
                                 RAM=5000, # MB
                                 run_ceiling_hours=24,
                                 run_ceiling_minutes=0,
                                 delim=',',v1=x,v2=y) # delimit records in the input files using commas because pickle uses newline characters.

                lsf_job.show_params()
                if loni_job:
                    if x == 0:
                        vals = {"jid":_job_group_id,"dep":_ids}
                        _last = _ids
                    else:
                        vals = {"jid":_job_group_id,"dep":_last}
                        _last = _ids
                    lsf_job.run(**vals)
                else:
                    # vals = {"dep":_trio[0]}
                    if x == 0:
                        vals = {"dep":_ids}
                        _last = _ids
                    else:
                        vals = {"dep":_last}
                        _last = _ids
                    lsf_job.run(**vals)

                x = y + 1
                if y + _node_count > len(trios) - 1:
                    y = len(trios) - 1 
                    x = y
                else:
                    y = y + _node_count
                _count = _count + 1
                if _count - 1 == 0:
                    triosUsedChecker.write("\n\nlen(trios): {0}, x: {1}, y: {2}, _loops: {3}, dep: {4}\noperands: {5}\n\n{6}".format(len(trios),0,_node_count,_loops,depp,str(lsf_job.operands),_ids[len(_ids)-1][0]))
                else:
                    triosUsedChecker.write("\n\nlen(trios): {0}, x: {1}, y: {2}, _loops: {3}, dep: {4}\noperands: {5}\n\n{6}".format(len(trios),x,y,_loops,depp,str(lsf_job.operands),_ids[len(_ids)-1][0]))

            triosUsedChecker.close()
    else:
        trios = [pickle.dumps(trio) for trio in trios]#moved here 5-25-17
        triosUsedChecker.write("\n\ninside else".format(len(trios)))
        lsf_job = LSF.LSFjob(worker_script='workerWGender.py',
                         operands=trios,
                         node_count=len(trios),
                         RAM=5000, # MB
                         run_ceiling_hours=24,
                         run_ceiling_minutes=0,
                         delim=',') # delimit records in the input files using commas because pickle uses newline characters.
        lsf_job.show_params()
        lsf_job.run()
