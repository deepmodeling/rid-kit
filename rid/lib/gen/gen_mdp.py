import re, os


def general_mdp(title, temperature=300, define=None, dt=0.002, nsteps=50000, frame_freq=500):
    if define is None:
        define = ''
    ret = ";  {}.mdp was generated.\n".format(title)
    # Preprocessor information: use cpp syntax. e.g.: -I/home/joe/doe -I/home/mary/roe
    ret += "include                  = \n"
    # e.g.: -DPOSRES -DFLEXIBLE (note these variable names are case sensitive)
    ret += "define                   = {} \n".format(define)
    ret += "integrator               = md\n"
    # Start time and timestep in ps
    ret += "tinit                    = 0\n"
    ret += "dt                       = {}\n".format(dt)
    ret += "nsteps = {}\n".format(nsteps)
    ret += "init-step                = 0\n"
    ret += "simulation-part          = 1 \n"
    ret += "comm-mode                = Linear\n"
    ret += "nstcomm                  = 1000\n"
    ret += "comm-grps                = \n"
    ret += "nstxout = {}\n".format(frame_freq)
    ret += "nstvout = {}\n".format(frame_freq)
    ret += "nstfout = {}\n".format(frame_freq)
    ret += "nstlog                   = 0\n"
    ret += "nstcalcenergy            = 10000\n"
    ret += "nstenergy = {}\n".format(frame_freq)
    ret += "nstxtcout = {}\n".format(frame_freq)
    ret += "xtc-precision            = 1000\n"
    ret += "xtc-grps                 = \n"
    ret += "energygrps               = \n"
    ret += "cutoff-scheme            = Verlet\n"
    ret += "nstlist                  = 10\n"
    ret += "ns-type                  = Grid\n"
    ret += "pbc                      = xyz\n"
    ret += "periodic-molecules       = no\n"
    ret += "verlet-buffer-drift      = -1\n"
    ret += "rlist                    = 1\n"
    ret += "coulombtype              = pme\n"
    ret += "rcoulomb                 = 0.9\n"
    ret += "epsilon-r                = 1\n"
    ret += "epsilon-rf               = 80\n"
    ret += "vdw-type                 = Cut-off\n"
    ret += "vdw-modifier             = Potential-shift-Verlet\n"
    ret += "rvdw                     = 0.9\n"
    ret += "DispCorr                 = EnerPres\n"
    ret += "table-extension          = 1\n"
    ret += "energygrp-table          = \n"
    ret += "fourierspacing           = 0.12\n"
    ret += "fourier-nx               = 0 \nfourier-ny               = 0 \nfourier-nz               = 0\n"
    ret += "pme-order                = 4\n"
    ret += "ewald-rtol               = 1e-05\n"
    ret += "ewald-geometry           = 3d\n"
    ret += "epsilon-surface          = 0\n"
    ret += "optimize-fft             = no\n"
    # Temperature
    ret += "tcoupl                   = v-rescale\n"
    ret += "tc-grps                  = water non-water\n"
    ret += "tau-t                    = 0.2 0.2\n"
    ret += "ref-t                    = {} {}\n".format(temperature, temperature)
    # Pressure
    ret += "pcoupl                   = parrinello-rahman\n"
    ret += "pcoupltype               = Isotropic\n"
    ret += "tau-p                    = 1.5\n"
    ret += "compressibility          = 4.5e-5\n"
    ret +="ref-p                    = 1.0\n"
    ret += "refcoord-scaling         = No\n"
    ret += "gen-vel                  = no\n"
    ret += "gen-temp                 = {}\n".format(temperature)
    ret += "gen-seed                 = 173529\n"
    ret += "constraints              = all-bonds\n"
    ret += "constraint-algorithm     = Lincs\n"
    ret += "continuation             = no\n"
    ret += "Shake-SOR                = no\n"
    ret += "shake-tol                = 0.0001\n"
    ret += "lincs-order              = 4\n"
    ret += "lincs-iter               = 1\n"
    ret += "lincs-warnangle          = 30\n"
    ret += "morse                    = no\n"
    ret += "energygrp-excl           = \n"
    ret += "nwall                    = 0\n"
    ret += "wall-type                = 9-3\n"
    ret += "wall-r-linpot            = -1\n"
    ret += "wall-atomtype            = \n"
    ret += "wall-density             = \n"
    ret += "wall-ewald-zfac          = 3\n"
    ret += "pull                     = no\n"
    ret += "rotation                 = no\n"
    ret += "user1-grps               =  \nuser2-grps               =  \nuserint1                 = 0 \nuserint2                 = 0 \nuserint3                 = 0 \nuserint4                 = 0\n"
    ret += "userreal1                = 0\n"
    ret += "userreal2                = 0\n"
    ret += "userreal3                = 0\n"
    ret += "userreal4                = 0\n"
    return ret


def gen_grompp_bias (nsteps, frame_freq, title='bias_md', temperature=300, define=None, dt=0.002) :
    # re.sub (pattern, subst, file_string)
    ret = general_mdp(title, temperature=temperature, nsteps=nsteps, frame_freq=frame_freq, dt=dt)
    return ret


def gen_grompp_res (nsteps, frame_freq, title='res_md', temperature=300, define="-DPOSRE", dt=0.002) :
    ret = general_mdp(title, temperature=temperature, define=define, nsteps=nsteps, frame_freq=frame_freq, dt=dt)
    re.sub("nstxout.*=.*", "nstxout = %d" % 0, ret)
    re.sub("nstvout.*=.*", "nstvout = %d" % 0, ret)
    re.sub("nstfout.*=.*", "nstfout = %d" % 0, ret)
    re.sub("nstenergy.*=.*", "nstenergy = %d" % 0, ret)
    return ret

def make_grompp(out_path, mdp_type, nsteps, frame_freq, title=None, temperature=300, define=None, dt=0.002):
    if mdp_type == 'bias':
        if os.path.basename(out_path) == '':
            out_path = os.path.abspath(out_path) + "/grompp.mdp"
        if title is None:
            title = 'bias_md'
        ret = gen_grompp_bias (nsteps, frame_freq, title=title, temperature=temperature, define=define, dt=dt)
    
    if mdp_type == 'res':
        if os.path.basename(out_path) == '':
            out_path = os.path.abspath(out_path) + "/grompp_restraint.mdp"
        if title is None:
            title = 'res_md'
        ret = gen_grompp_res (nsteps, frame_freq, title=title, temperature=temperature, define=define, dt=dt)
    
    with open(out_path, 'w') as mdp:
        mdp.write(ret)


def modify_grompp_bias (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % frame_freq)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % frame_freq)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % frame_freq)
    replace (gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)    

def modify_grompp_res (gro_file, nsteps, frame_freq) :
    replace (gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace (gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace (gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace (gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace (gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace (gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)
    

if __name__ == '__main__':
    con = gen_mdp('mdout', 320)
    with open('./test.mdp', 'w') as mdp:
        mdp.write(con)