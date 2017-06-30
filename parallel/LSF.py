import os
import sys

"""
    worker script requirements:
    --------------------------------------------------
    The worker script has 2 requirements:
      1) It must accept the path to the input file as command line argument.
      2) All output must be through stdout (just use a python print statement).
         The LSF process will redirect this to a useful output file.

    i.e. a call to the worker script should look like this:
    $ python worker_script.py input.txt
"""

class LSFjob:
    def __init__(self, worker_script, operands, node_count, RAM, run_ceiling_hours=0, run_ceiling_minutes=0, delim='\n',v1=0,v2=0):
        # a list of items to be operated on
        self.operands = operands
        # the hours componant of the enforced max runtime
        self.run_ceiling_hours = run_ceiling_hours
        # the minutes componant of the enforced max runtime
        self.run_ceiling_minutes = run_ceiling_minutes
        # the RAM (MB) allocated to each node
        self.RAM = RAM
        # the number of nodes
        self.node_count = node_count
        # the path to the worker script
        self.worker_script = worker_script
        # a file delimiter for the records in the input file
        self.delim = delim
        # a directory to hold input and output files
        self.input_output_dir = './lsf_data'

        self._check_params() # check the input types and bounds
        #######################################################################################
        #5-26-17
        self.startx = v1; self.endy=v2; self.bulk_job = False; self.bulk_file_numbers = None;
        #######################################################################################
        self._populate_data_directory()


    def _check_params(self):
        # check the input parameter types and bounds
        if isinstance(self.run_ceiling_hours, (int, float, long)) and self.run_ceiling_minutes >= 0:
            self.run_ceiling_hours = int(self.run_ceiling_hours)
        else:
            raise Exception("(input error) run_ceiling_hours must be a non-negative integer.")

        if isinstance(self.run_ceiling_minutes, (int, float, long)) and self.run_ceiling_minutes >= 0:
            self.run_ceiling_minutes = int(self.run_ceiling_minutes)
            if self.run_ceiling_minutes > 60:
                self.run_ceiling_minutes = 60
        else:
            raise Exception("(input error) run_ceiling_minutes must be a non-negative integer.")

        if isinstance(self.RAM, (int,float,long)) and self.RAM >= 0:
            self.RAM = int(self.RAM)
        else:
            raise Exception("(input error) RAM must be a non-negative integer.")

        if isinstance(self.node_count, (int,float,long)) and self.node_count > 0 and self.node_count <= 70:
            self.node_count = int(self.node_count)
        else:
            raise Exception("(input error) node_count must be an integer between 1 and 70, inclusive.")

        if not os.path.isfile(self.worker_script):
            raise Exception("(input error) Can't find the file, \'%s\'" % (self.worker_script))

        if len(self.operands) < self.node_count:
            raise Exception("The number of operands must be >= to the number of nodes.")

    def _calculate_loads(self):
        # create a balanced list of the number of operands each node will receive.
        operand_count = len(self.operands)
        base_load_size = operand_count/self.node_count
        large_load_count = operand_count%self.node_count
        larger_loads = [base_load_size+1 for i in range(large_load_count)]
        smaller_loads = [base_load_size for i in range(self.node_count - large_load_count)]
        return larger_loads + smaller_loads

    def _input_path(self,i):
        return os.path.join(self.input_output_dir,"input_%i.dat" % (i))

    def _output_path(self,i):
        return os.path.join(self.input_output_dir,"output_%i.dat" % (i))

    def _populate_data_directory(self):
        # if the I/O directory doesn't exists, create it.
        if not os.path.isdir(self.input_output_dir):
            os.mkdir(self.input_output_dir)
        # get a list of load sizes for each node
        load_sizes = self._calculate_loads()
        # create a copy of the operands list so that we can pop off consumed elements
        operands = self.operands[:]
        ##################################################################################
        #5-26-17
        if self.startx != self.endy:
            self.bulk_job = True
            if self.startx == 0:
                self.bulk_file_numbers = range(self.startx,self.endy)
            else:
                self.bulk_file_numbers = range(self.startx-1,self.endy)
        ##################################################################################
        for i in range(self.node_count):
            # get the load size for this node
            load_size = load_sizes[i]
            # extract a list of operands for this node based on the prescribed load size
            load = operands[:load_size]
            # "pop" the consumed operands off the big list
            del operands[:load_size]
            # create and populate files
        ##################################################################################
        #5-25-17
            if self.bulk_job:
                with open(self._output_path(self.bulk_file_numbers[i]),'w') as f:
                    pass
                with open(self._input_path(self.bulk_file_numbers[i]),'w') as f:
                    f.write(self.delim.join([str(operand) for operand in load]))
        ##################################################################################
            else:
                with open(self._output_path(i),'w') as f:
                    pass
                with open(self._input_path(i),'w') as f:
                    f.write(self.delim.join([str(operand) for operand in load]))

    def show_params(self):
        print "Params Summary"
        print "--------------------------------------------"
        print "Run time ceiling:   %i:%i (hours:minutes)" % (self.run_ceiling_hours,self.run_ceiling_minutes)
        print "RAM per node:       %i MB" % (self.RAM)
        print "Number of nodes:    %i" % (self.node_count)
        print "Worker script path: %s" % (self.worker_script)
        print "Number of operands: %i" % (len(self.operands))
        print "--------------------------------------------"
        if self.startx == 0:
            print "range: %i-%i"%(self.startx,self.endy)
        else:
            print "range: %i-%i"%(self.startx-1,self.endy)
    
    def run(self,**kwargs):
        ###############
        #5-24-17
        ###############
        if kwargs:
            if self.bulk_job:
                for c,i in enumerate(self.bulk_file_numbers):
                    # if i == 0:#testing to make sure the next job starts after the first jobs finishes 5-30-17 
                    if i < self.node_count:
                        if "jid" in kwargs.keys() and kwargs['jid'] is not None:
                            add = {"J": str(kwargs['jid']) + "-" + str(i)}
                            cmd = self.lsf_command_wrapper(i,**add)
                        else:
                            add = {"J": "nonloni-" + str(i)}
                            cmd = self.lsf_command_wrapper(i,**add)
                        os.popen(cmd)
                        print cmd

                    elif "dep" in kwargs.keys() and kwargs['dep'] is not None:
                        try:
                            # _job_done_file = "{0}_output/{0}-trio-validation-complete.txt".format(kwargs['dep'][c])#works
                            _trio_id = "{0}".format(kwargs['dep'][c][0][0:7])
                        except Exception as _bad_dep_index:
                            print("Exception: {0}\nDep: {1}".format(_bad_retry_dep_index,kwargs['dep']))
                            sys.exit("Bad dep index! Safe exit...")
                        if "jid" in kwargs.keys() and kwargs['jid'] is not None:
                            # add = {"df":_trio_id,"J": str(kwargs['jid']) + "-" + str(i), "w": str(kwargs['jid']) + "-" + str(i-(self.node_count))}
                            add = {"df":_trio_id,"J": str(kwargs['jid']) + "-" + str(i)}
                            if len(kwargs['dep'][c]) > 3:
                                add.update({"multi-key-dep":kwargs['dep'][c]})
                            cmd = self.lsf_command_wrapper(i,**add)
                        else:
                            # add = {"df":_trio_id,"J": "nonloni-" + str(i), "w": "nonloni-"+ str(i-(self.node_count))}
                            add = {"df":_trio_id,"J": "nonloni-" + str(i)}
                            if len(kwargs['dep'][c]) > 3:
                                add.update({"multi-key-dep":kwargs['dep'][c]})
                            cmd = self.lsf_command_wrapper(i,**add)
                        os.popen(cmd)
                        print "dep: {0}\ntrio dependency: {1}\nlength trio dep: {2}".format(kwargs['dep'],_trio_id,len(kwargs['dep'][c]))
                        print cmd
        else:
            for i in range(self.node_count):
                cmd = self.lsf_command_wrapper(i)
                os.popen(cmd)
                print cmd

    def lsf_command_wrapper(self,c,**extraOptions):
        cmd = "bsub -W {0}:{1} -M {2} -o {3} \"python {4} {5} {6}\"".format(self.run_ceiling_hours,
                                                                   self.run_ceiling_minutes,
                                                                   self.RAM,
                                                                   self._output_path(c),
                                                                   self.worker_script,
                                                                   self._input_path(c),
                                                                   len(self.operands))

        if extraOptions:
            if "J" in extraOptions.keys():
                cmd += " -J \"{0}\"".format(extraOptions["J"])
            if "df" in extraOptions.keys():
                if "multi-key-dep" not in extraOptions.keys():
                    cmd = "while ! (test -e \"{0}_output/{0}-done.txt\" && test -e \"{0}_output/{0}-01-done.txt\" && test -e \"{0}_output/{0}-02-done.txt\"); do sleep 90; done;\n{1}".format(extraOptions['df'],cmd)
                else:
                    cmd = "while ! ("
                    for i,df in enumerate(extraOptions["multi-key-dep"]):
                        if i < len(extraOptions["multi-key-dep"]) - 1:
                            cmd += "test -e \"{0}_output/{1}-{2}-done.txt\" && ".format(extraOptions['df'],df,i)
                        elif i == len(extraOptions["multi-key-dep"]) - 1:
                            cmd += "test -e \"{0}_output/{1}-{2}-done.txt\"); do sleep 90; done;\n".format(extraOptions['df'],df,i)

            else:
                cmd += "; " 
        cmd += "sleep 10"
        return cmd
