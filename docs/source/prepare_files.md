# Preparing files for input to `RiD-kit`
Currently Rid-kit manage files by their extensions, so you have to provide files with the correct extensions in the input directory. These includes the following:
## For `Gromacs` simulation

`.gro files` represent the initial configurations used by Gromacs run.

`.top file` represent the topology information used by Gromacs run.

`.ff directory (optional)` represents the forcefield information used by Gromacs, all the forcefield and restraint information not implemented in Gromacs natively should be included in this directory. Otherwise you can manually specify the forcefield information in the topology file.

`.ndx file (optional)` represents the index file used by Gromacs, in situations where we want the control of different parts of the system.

`.pb files (optional)` represents either the free energy model produced by Rid run, or the Deep Potential model used to run Gromacs. Specify different models in `rid.json` according to [Rid configuration](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_configuration.md).

`CV file (any extension or no extension) and its related .pdb files (optional)`: The CV file is used to customize CV definition, other files related to CV definition should has `pdb` extensions. The rules for defining customized CV can be found in [Rid configuration](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_configuration.md).

`.npy file (optional)` represents the CV-mean forces data used to train free energy model, this file is used when you want to use data from a previous RiD run.

## For `Lammps` simulation

`.lmp files` represent the initial configurations used by Lammps run.

`input file (any extension or no extension)` represents the input file used by Lammps run.

`.pb files (optional)` represents either the free energy model produced by Rid run, or the Deep Potential model used to run Gromacs. Specify different models in `rid.json` according to [Rid configuration](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_configuration.md).

`CV file (any extension or no extension) and its related .pdb files (optional)`: The CV file is used to customize CV definition, other files related to CV definition should has `pdb` extensions. The rules for defining customized CV can be found in [Rid configuration](https://github.com/PKUfjh/rid-kit/blob/dflow/docs/source/rid_configuration.md).

`.npy file (optional)` represents the CV-mean forces data used to train free energy model, this file is used when you want to use data from a previous RiD run.