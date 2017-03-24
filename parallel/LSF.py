import os

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
    def __init__(self, worker_script, operands, node_count, RAM, run_ceiling_hours=0, run_ceiling_minutes=0, delim='\n'):
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
        for i in range(self.node_count):
            # get the load size for this node
            load_size = load_sizes[i]
            # extract a list of operands for this node based on the prescribed load size
            load = operands[:load_size]
            # "pop" the consumed operands off the big list
            del operands[:load_size]
            # create and populate files
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

    def run(self):
        for i in range(self.node_count):
            cmd = 'bsub -W %i:%i -M %i -o %s "python %s %s %s"' %(self.run_ceiling_hours,
                                                               self.run_ceiling_minutes,
                                                               self.RAM,
                                                               self._output_path(i),
                                                               self.worker_script,
                                                               self._input_path(i),
                                                               len(self.operands))
            os.popen(cmd)
            print cmd
