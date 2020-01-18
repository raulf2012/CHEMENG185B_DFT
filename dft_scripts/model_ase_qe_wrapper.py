#!/usr/bin/env python

"""Simple DFT optimization example.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

from ase import io

from ase.calculators.espresso import Espresso
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

calc = Espresso(
    kpts=(3, 3, 3),
    ecutwfc=400 * eV2Ry,

    # 'smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'from_input',
    occupations="smearing",
    # 'gaussian', 'methfessel-paxton', 'marzari-vanderbilt', 'fermi-dirac'
    smearing="gaussian",
    degauss=0.05 * eV2Ry,

    tstress=True,
    tprnfor=True,

    pseudopotentials=pseudopotentials,
    )
atoms.set_calculator(calc)
# __|


atoms.get_potential_energy()

# try:
#     atoms.get_forces()
#     print("Getting forces worked!!")
# except:
#     print("getting forces didn't work")

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
# cubic lattic constant
# print((8*rocksalt.get_volume()/len(rocksalt))**(1.0/3.0))

# __|
