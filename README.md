# Protein-ligand Complex Parameterization 

This program is to automatically generate parameterized protein-ligand complex for Molecular Dynamics simulation using antechamber and gromacs with Amber forcefield.  


## How to run the code. 

### Inputs
The code assumes the input files are seperated into protein file and ligands poses, both in PDB formats, which are obtained from WF0 outputs. It will require a script from WF0 to convert their outputs to PDB formate. The protein PDB is specified at line 20 and ligands at line 16. The ligand PDBs are fetched through `glob` that works similarly as `ls`.  

### Running the code 
The coda is simply running 
```
source activate entk_py3
export RADICAL_PILOT_DBURL=mongodb://hyperrct:v9gpU8gU5Wxz2C88@129.114.17.185:27017/hyperrct
export RMQ_HOSTNAME=129.114.17.185
export RMQ_PORT=5672
export RMQ_USERID=hyperrct
export RMQ_PASSWD=h1p3rrc7
export RADICAL_PROFILE=True
export RADICAL_PILOT_PROFILE=True
export RADICAL_ENTK_PROFILE=True
python summit_param.py
```
after all the input path is specified in the script. 


### Outputs 
The output will be in N created folders, where N is number of ligands plus 1 protein. It should be ready for WF2 in the version running adrp. 

## Dependencies for ff setup 
This is the most tricky part as we are missing some of them on Summit. 
1. ambertools 
1. openbabel
2. MDAnalysis
2. OpenMM
3. tqdm
4. subprocess
5. parmed
6. numpy
7. h5py
8. tempfile 

Most of these should trivial except for ambertools and openbabel. I have never tried to install them on Summit, nor there is compiled version from conda. It was normally run on my laptop or other Linux machines.  
So the python is 3.6 and most dependencies should be available here. (/ccs/home/hm0/.conda/envs/py3) Since both ambertools and openbabel can be installed independent of conda environment, `py` env can be loaded for all the python dependencies. 
The `environment.yml` can be a good reference if the list is missing any package. 

### Installations 
1. ambertools 
    PowerPC version of ambertools was available from the `bioconda` channel. 
    ```
        conda install -c bioconda ambertools
    ```
2. openbabel 
    OpenBabel was compiled from source with `cmake`. 
    ```
        cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.conda/envs/py3/ -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
    ```
