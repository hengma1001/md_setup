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
It requires pdb file as input format. The amber function can automatically parameterize small molecule ligand with antechamber. The gmx function needs an additional entry to specify the force field path. 

### Running the code 
The coda is simply running 
```
python param_amber.py
```
after all the input path is specified in the script. 

