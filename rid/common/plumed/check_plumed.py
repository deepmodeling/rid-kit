from typing import Dict


def check_deepfe_input(inputs: Dict):
    assert "mode" in inputs.keys(), "Please specify a mode for making plumed files."
    if mode == "torsion":
        assert "conf" in inputs.keys(), "Please provide a conformation to make torsions."
        assert "selected_id" in inputs.keys(), "Please provide selected residues"
    elif mode == "custom":
        assert "cv_file" in inputs.keys()
    else:
        raise RuntimeError("Unknown mode to make plumed files.")
    assert float(inputs["trust_lvl_1"]) < float(inputs["trust_lvl_2"])
    assert "stride" in inputs.keys()
    assert "output" in inputs.keys()
