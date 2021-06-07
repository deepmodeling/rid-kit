from rid.lib.utils import replace


def make_grompp_bias(gro_file, nsteps, frame_freq):
    replace(gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace(gro_file, "nstxout.*=.*", "nstxout = %d" % frame_freq)
    replace(gro_file, "nstvout.*=.*", "nstvout = %d" % frame_freq)
    replace(gro_file, "nstfout.*=.*", "nstfout = %d" % frame_freq)
    replace(gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace(gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)


def make_grompp_res(gro_file, nsteps, frame_freq):
    replace(gro_file, "nsteps.*=.*", "nsteps = %d" % nsteps)
    replace(gro_file, "nstxout.*=.*", "nstxout = %d" % 0)
    replace(gro_file, "nstvout.*=.*", "nstvout = %d" % 0)
    replace(gro_file, "nstfout.*=.*", "nstfout = %d" % 0)
    replace(gro_file, "nstxtcout.*=.*", "nstxtcout = %d" % frame_freq)
    replace(gro_file, "nstenergy.*=.*", "nstenergy = %d" % frame_freq)
