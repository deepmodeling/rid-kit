dihedral_name = "dih-{resid:03d}-{angid:02d}"
dihedral_def_from_atoms = "{name}: TORSION ATOMS={a1},{a2},{a3},{a4}"
distance_name = "dis-{atomid1:05d}-{atomid2:05d}"
distance_def_from_atoms = "{name}: DISTANCE ATOMS={a1},{a2}"
deepfe_def = "dpfe: DEEPFE TRUST_LVL_1={trust_lvl_1} TRUST_LVL_2={trust_lvl_2} MODEL={model} ARG={arg}"
print_def = "PRINT STRIDE={stride} ARG={arg} FILE={file}"
restraint_def = "{name}: RESTRAINT ARG={arg} KAPPA={kappa} AT={at}"
restraint_prefix = "res"
upper_def = "UPPER_WALLS ARG={arg} AT={at} KAPPA={kappa} LABEL={name}"
lower_def = "LOWER_WALLS ARG={arg} AT={at} KAPPA={kappa} LABEL={name}"