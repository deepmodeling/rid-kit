# rid task name format
walker_tag_fmt = "{idx:03d}"
init_conf_name = "conf_{idx:03d}.gro"
explore_task_pattern = "{:03d}"
explore_task_file = "explore_{walker:03d}.pkl"
culster_selection_data_name = "cls_sel.out"
culster_selection_index_name = "cls_sel.ndx"
sel_ndx_name = "sel.ndx"
cv_init_label = "cv_init_{idx:d}.out"
model_devi_name = "model_devi.txt"
label_task_pattern = "{:03d}"
force_out = "forces.out"
data_new = "data.new"
data_old = "data.old"
data_raw = "data.raw"
block_tag_fmt = "iter-{idx_iter:03d}"

# PLUMED2 file names
plumed_input_name = "plumed.dat"
plumed_bf_input_name = "plumed_bf.dat"
plumed_restraint_input_name = "plumed_restraint.dat"
plumed_output_name = "plm.out"

# Gromacs file names
gmx_conf_name = "conf.gro"
gmx_top_name = "topol.top"
gmx_mdp_name = "grompp.mdp"
gmx_tpr_name = "topol.tpr"
gmx_grompp_log = "gmx_grompp.log"
gmx_mdrun_log = "md.log"
restraint_md_mdp_name = "grompp_restraint.mdp"
trr_name = "traj.trr"
xtc_name = "traj_comp.xtc"
sel_gro_name = "conf_{idx:d}.gro"
sel_gro_name_gmx = "conf_.gro"
gmx_conf_out = "confout.gro"


# Tensorflow files
tf_model_name = "model_{tag}.pb"
model_tag_fmt = "{idx:03d}"
N_grid = 100


# Units
kbT = (8.617343E-5) * 300
beta = 1.0 / kbT
f_cvt = 96.485
inverse_f_cvt = 1 / f_cvt

# precision
cls_ndx_precision = "%d"
cls_sel_precision = "%.6e"
model_devi_precision = "%.6e"