# Protein-ligand Complex Parameterization 

This program is to automatically generate parameterized protein-ligand complex for Molecular Dynamics simulation using antechamber or gromacs `pdb2gmx`.

## Before running the code

The gromacs portion of the code doesn't work on Mac for `echo` behaves differently in Mac from Linux. The `pdb2gmx` manual entry fails to register in the Mac environment. 

### Set up the env

The environment can be easily built with `conda` via

```bash
conda env create -f env.yml
```

### Load the env 

The conda environment should be available to load with 

```bash
source activate MD_ff
```



## Run the parameterization

The code takes in pdb files and builds topology according to the chosen force field. Examples can be found [here][https://github.com/hengma1001/complex_sim/tree/master/examples]. 

### Inputs
The code assumes the input files are seperated into protein file and ligands poses, both in PDB formats, which are obtained from WF0 outputs. It will require a script from WF0 to convert their outputs to PDB formate. The protein PDB is specified at line 20 and ligands at line 16. The ligand PDBs are fetched through `glob` that works similarly as `ls`.  

### Running the code 
The coda is simply running 
```
python parameterize.py
```
after all the input path is specified in the script. 


### Outputs 
The output will be in N created folders, where N is number of ligands plus 1 protein. It should be ready for WF2 in the version running adrp. 

