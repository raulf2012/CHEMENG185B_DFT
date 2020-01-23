#!/usr/bin/env python

"""Simple DFT optimization example.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

from ase import io

from espresso import espresso

# from ase.calculators.espresso import Espresso
from ase.optimize import LBFGS

# from ase.build import bulk
# from ase.constraints import UnitCellFilter
#__|


# Read atoms object
atoms = io.read("init.traj")
eV2Ry = 0.073498618


#| - QE Calculator
pseudopotentials = {
    "C": "C.UPF",
    "Fe": "Fe.UPF",
    "H": "H.UPF",
    }

# calc = Espresso(
# calc = espresso(
#     kpts=(3, 3, 3),
#     ecutwfc=400 * eV2Ry,
#
#     # 'smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'from_input',
#     occupations="smearing",
#     # 'gaussian', 'methfessel-paxton', 'marzari-vanderbilt', 'fermi-dirac'
#     smearing="gaussian",
#     degauss=0.05 * eV2Ry,
#
#     tstress=True,
#     tprnfor=True,
#
#     pseudopotentials=pseudopotentials,
#     )


calc = espresso(
    pw=300,
    dw=3000,
    xc='PBE',
    kpts=[1, 1, 1],
    sigma=0.1,
    smearing='mv',
    spinpol=False,
    convergence={
        'energy': 0.0001,
        'mixing': 0.1,
        'nmix': 10,
        'maxsteps': 300,
        'diag': 'david',
        },

    output={
        'avoidio': True,
        'removewf': True,
        'wf_collect': False,
        'removesave': True,
        },
    outdir='calcdir',
    )

atoms.set_calculator(calc)
# __|


atoms.get_potential_energy()


#| - QuasiNewton
# qn = QuasiNewton(
#     atoms,
#     # trajectory="out_opt.traj",
#     logfile="qn.log",
#     )
#
# if traj is not None:
#     qn.attach(traj)  # COMBAK Test feature (restarting traj files)
#
# qn.run(
#     fmax=0.005,
#     steps=100,
#     )
# __|


opt = LBFGS(
    atoms,
    restart="restart.traj",
    logfile="log",
    trajectory="out.traj",
    maxstep=None,
    memory=100,
    damping=1.0,
    alpha=70.0,
    use_line_search=False,
    master=None,
    force_consistent=None,
    )
opt.run(fmax=0.005)

atoms.write("final.traj")



#| - __old__
# try:
#     atoms.get_forces()
#     print("Getting forces worked!!")
# except:
#     print("getting forces didn't work")


# cubic lattic constant
# print((8*rocksalt.get_volume()/len(rocksalt))**(1.0/3.0))
# __|
