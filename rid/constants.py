# rid task name format
argo_namespace = "argo"
walker_tag_fmt = "{idx:03d}"
init_conf_gmx_name = "conf_{idx:03d}.gro"
init_conf_lmp_name = "conf_{idx:03d}.lmp"
init_input_name = "input_{idx:03d}.lammps"
explore_task_pattern = "{:03d}"
explore_task_file = "explore_{walker:03d}.pkl"
cluster_selection_data_name = "cls_sel.out.npy"
cluster_selection_index_name = "cls_sel.ndx.npy"
sel_ndx_name = "sel.ndx.npy"
cv_init_label = "cv_init_{walker:03d}_{idx:d}.out"
model_devi_name = "model_devi.txt"
label_task_pattern = "{:03d}"
cv_force_out = "cv_forces.out"
data_new = "data.new.npy"
data_old = "data.old.npy"
data_raw = "data.raw.npy"
block_tag_fmt = "iter-{idx_iter:03d}"

# PLUMED2 file names
plumed_input_name = "plumed.dat"
plumed_bf_input_name = "plumed_bf.dat"
plumed_restraint_input_name = "plumed_restraint.dat"
plumed_output_name = "plm.out"
center_out_name = "centers.out"

# Gromacs file names
gmx_conf_name = "conf.gro"
gmx_top_name = "topol.top"
gmx_idx_name = "index.ndx"
gmx_mdp_name = "grompp.mdp"
gmx_tpr_name = "topol.tpr"
gmx_grompp_log = "gmx_grompp.log"
gmx_mdrun_log = "md.log"
restraint_md_mdp_name = "grompp_restraint.mdp"
gmx_trr_name = "traj.trr"
gmx_xtc_name = "traj_comp.xtc"
sel_gro_name = "conf_{walker:03d}_{idx:d}.gro"
sel_gro_name_gmx = "conf_.gro"
gmx_conf_out = "confout.gro"
gmx_coord_name = "coord.xvg"
gmx_force_name = "force.xvg"

# Lammps file names
lmp_conf_name = "conf.lmp"
lmp_input_name = "input.lammps"
lmp_mdrun_log = "log.lammps"
sel_lmp_name = "conf_{walker:03d}_{idx:d}.lmp"
lmp_conf_out = "confout.lmp"


# Tensorflow files
tf_model_name = "model_{tag}.pb"
train_log = "log_{tag}"
dp_model_name = "dp.pb"
model_tag_fmt = "{idx:03d}"
N_grid = 100

# Dp config file
dp_config_name = "dp_config"
dp_sel_ndx = "dp_sel.ndx"

# figure file name
mf_fig = "mf_average.png"
bias_fig = "bias.png"
model_devi_fig = "model_devi.png"
new_model_devi_fig = "new_model_devi.png"
cluster_fig = "cluster.png"
train_fig = "train_{tag}.png"
dp_model_devi_fig = "dp_model_devi.png"
mf_std_fig = "mf_std.png"

# MCMC file name
mcmc_1cv_dir_name = "mcmc_1cv"
mcmc_2cv_name = "mcmc_2cv.dat"
mcmc_1cv_name = "mcmc_1cv_{tag}.dat"
mcmc_2cv_fig = "mcmc_2cv.png"
mcmc_1cv_fig = "mcmc_1cv_{tag}.png"
mcmc_2cv_fig_separate = "mcmc_2cv_{tag}.png"

# Units
kb = 8.617333E-5 # ev
kbT = (8.617333E-5) * 300
beta = 1.0 / kbT
f_cvt = 96.485 # ev to kj/mol
inverse_f_cvt = 1 / f_cvt
kcal2kj = 4.184

# precision
model_devi_precision = "%.6e"