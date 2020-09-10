import os, json, time 
import glob
from radical.entk import Pipeline, Stage, Task, AppManager

# ------------------------------------------------------------------------------
# Set default verbosity

if os.environ.get('RADICAL_ENTK_VERBOSE') is None:
    os.environ['RADICAL_ENTK_REPORT'] = 'True'

# Assumptions:
# - # of MD steps: 2
# - Each MD step runtime: 15 minutes
# - Summit's scheduling policy [1]
#
# Resource rquest:
# - 4 <= nodes with 2h walltime.
#
# Workflow [2]
#
# [1] https://www.olcf.ornl.gov/for-users/system-user-guides/summit/summit-user-guide/scheduling-policy
# [2] https://docs.google.com/document/d/1XFgg4rlh7Y2nckH0fkiZTxfauadZn_zSn3sh51kNyKE/
#
'''
export RADICAL_PILOT_DBURL=mongodb://hyperrct:v9gpU8gU5Wxz2C88@129.114.17.185:27017/hyperrct
export RMQ_HOSTNAME=129.114.17.185
export RMQ_PORT=5672
export RMQ_USERID=hyperrct
export RMQ_PASSWD=h1p3rrc7
export RADICAL_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export RADICAL_ENTK_PROFILE=True
'''
#

base_path = os.path.abspath('.') 
conda_path = '/ccs/home/hm0/.conda/envs/py3' 

prot_file = os.path.join(base_path, 'pdbs/prot.pdb') 
lig_files = glob.glob(
        os.path.join(base_path, 'pdbs/input_*.pdb')
        )

N_jobs = 1 + len(lig_files) 
print(prot_file, lig_files, N_jobs) 

hrs_wt = 6 
queue = 'batch'
proj_id = 'med110'

batch_size = 256

class Parameterization:
    """
    """

    def __init__(
            self, 
            prot_file, 
            lig_files,
            ): 
        """
        Parameterization of protein and compounds 

        Parameters
        ----------
        prot_file : str
            Path to the protein pdb
        lig_files : list
            List of ligand pdb files 
        """

        self.prot_file = prot_file
        self.lig_files = lig_files
        self.n_jobs = 1 + len(lig_files) 


    def param_pipeline(self): 
        """ 
        Function to parameterize all parameterization pipeline
        """ 
        p = Pipeline() 
        p.name = 'param' 
        s = Stage()
        s.name = 'param'

        # Aggregation task
        for i in range(self.n_jobs):
            t = Task()
            # https://github.com/radical-collaboration/hyperspace/blob/MD/microscope/experiments/MD_to_CVAE/MD_to_CVAE.py
            t.pre_exec = [] 
            t.pre_exec += ['. /sw/summit/python/2.7/anaconda2/5.3.0/etc/profile.d/conda.sh']
            t.pre_exec += ['conda activate %s' % conda_path] 
            t.pre_exec += ['cd %s' % base_path]
            t.executable = ['%s/bin/python' % conda_path]  # MD_to_CVAE.py
            t.arguments = [
                    '%s/parameterize.py' % base_path, 
                    '--prot', self.prot_file
                    ]
            if i > 0: 
                lig_file = self.lig_files[i-1]
                t.arguments += ['--lig', lig_file]

            # assign hardware the task 
            t.cpu_reqs = {
                    'processes': 1,
                    'process_type': None,
                    'threads_per_process': 2,
                    'thread_type': 'OpenMP'
                    }
            # Add the aggregation task to the aggreagating stage
            s.add_tasks(t)

        p.add_stages(s)
        return p


if __name__ == '__main__':

    # Create a dictionary to describe four mandatory keys:
    # resource, walltime, cores and project
    # resource is 'local.localhost' to execute locally
    res_dict = {
            'resource': 'ornl.summit',
            'queue'   : queue,
            'schema'  : 'local',
            'walltime': 60 * hrs_wt,
            'cpus'    : N_jobs * 2,
            # 'gpus'    : n_gpus,#6*2 ,
            'project' : proj_id
    }

    Param = Parameterization(prot_file, lig_files)
    p1 = Param.param_pipeline()
    # Create Application Manager
    # appman = AppManager()
    appman = AppManager(
            hostname=os.environ.get('RMQ_HOSTNAME'), 
            port=int(os.environ.get('RMQ_PORT')), 
            username=os.environ.get('RMQ_USERID'), 
            password=os.environ.get('RMQ_PASSWD'))
    appman.resource_desc = res_dict


    pipelines = [p1] #[]
    # pipelines.append(p1)
    # pipelines.append(p2)

    # Assign the workflow as a list of Pipelines to the Application Manager. In
    # this way, all the pipelines in the list will execute concurrently.
    appman.workflow = pipelines

    # Run the Application Manager
    appman.run()
