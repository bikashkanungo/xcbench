import sys
import math
from pyscf import gto, scf
from pyscf import dft
import numpy as np

atomSym = sys.argv[1]
basis = sys.argv[2]
lmap = {0: 'S', 1: 'P', 2: 'D', 3: 'F', 4: 'G', 5: 'H', 6: 'I'}
mol = gto.M(
    atom = f'''
    {atomSym}  0.   0.       0.
    ''',
    basis = basis,
    spin=1)
nbasis = mol.nao_nr()
print('nbasis', nbasis)
outF = open(atomSym + "_" + basis, "w")
print(atomSym, 1.0, file = outF)

nshells = mol.atom_nshells(0)
for i in range(nshells):
    l = mol.bas_angular(i)
    c = mol.bas_ctr_coeff(i)
    e = mol.bas_exp(i)
    nc = e.shape[0]
    nd = c.shape[1]
    lsym = lmap[l]
    for j in range(nd):
        print(lsym, nc, file = outF)
        for k in range(nc):
            print(e[k], c[k,j], file = outF)
        
outF.close()
