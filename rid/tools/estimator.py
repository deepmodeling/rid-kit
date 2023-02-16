import numpy as np

def pseudo_inv(r_cv, cv_init, selected_atomid_simple):
    A = np.zeros((len(r_cv),len(cv_init)))
    for index in range(len(cv_init)):
        atom_id1 = selected_atomid_simple[index][0]
        atom_id2 = selected_atomid_simple[index][1]
        A[atom_id1*3,index] = (r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        A[atom_id1*3+1,index] = (r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        A[atom_id1*3+2,index] = (r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
        A[atom_id2*3,index] = -(r_cv[atom_id1*3] - r_cv[atom_id2*3])/cv_init[index]
        A[atom_id2*3+1,index] = -(r_cv[atom_id1*3+1] - r_cv[atom_id2*3+1])/cv_init[index]
        A[atom_id2*3+2,index] = -(r_cv[atom_id1*3+2] - r_cv[atom_id2*3+2])/cv_init[index]
    
    U, S, Vh = np.linalg.svd(A, full_matrices = False)
    B = np.matmul(np.matmul(np.transpose(Vh),np.linalg.inv(np.diag(S))),np.transpose(U))
    
    return B