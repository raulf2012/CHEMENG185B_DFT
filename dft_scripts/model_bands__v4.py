#!/usr/bin/env python
#| - SLURM HEADER
#################
#SBATCH --job-name=test00
#################
#SBATCH --time=460:00
#################
#SBATCH --ntasks=8
#################
#SBATCH --nodes=1
#################
#SBATCH --mem-per-cpu=6G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END,FAIL
#################
#who to send email to; please change to your email
# #SBATCH  --mail-user=flores12@stanford.edu
#################
#a file for job output, you can check job progress
#SBATCH --output=job.out
#################
# a file for errors from the job
#SBATCH --error=job.err
#################
#__|


#| - Import Modules
import os
import sys

import shutil
import pickle

from ase import io
from ase.dft.kpoints import ibz_points, get_bandpath

from vossjo_espresso import espresso as Espresso
#__|

#| - Read atoms object
try:
    atoms = io.read("init.traj")
except:
    atoms = io.read("init.cif")
#__|

#| - Espresso calculator
regular_params = dict(
    pw=500,  # IMP | (ev)
    nbands=-50,
    kpts=(6, 6, 6),
    mode="scf",
    opt_algorithm="bfgs",
    convergence={
        'energy': 1e-4,
        'mixing': 0.1,
        'maxsteps': 300,
        'diag': 'david',
        },
    )

easy_params = dict(
    pw=50,  # IMP | (ev)
    nbands=-5,
    kpts=(1, 1, 1),
    mode="scf",
    opt_algorithm="bfgs",
    # calculation='scf',
    convergence={
        'energy': 1e-2,
        'mixing': 0.1,
        'maxsteps': 300,
        'diag': 'david',
        },
    )


# calc = espresso(
calc = Espresso(
    kptshift=(0, 0, 0),

    U={
        "V": 2.,
        #  "V": 2.5,
        #  "V": 3.,
        #  "V": 3.5,
        #  "V": 4.,
        #  "V": 4.5,

        "Bi": 0.,
        "O": 0.,
        },

    # ion_dynamics='bfgs',
    fmax=0.05,
    xc='PBE',  # IMP
    outdir="calcdir",
    U_projection_type='atomic',
    # site=site_config,
    tot_magnetization=-1,
    occupations='smearing',
    dipole={'status': False},
    field={'status': False},
    output={
        "avoidio": False,
        "removewf": False,
        "removesave": False,  # Changed this to true so it wouldn't make big files anymore
        "wf_collect": False,

        # 'disk_io': 'default',
        # # 'avoidio': True,
        # 'removewf': True,
        # 'removesave': True,
        # 'wf_collect': False,
        },
    smearing='fd',
    sigma=0.1,
    # #########################################################################
    alwayscreatenewarrayforforces=True,
    verbose='low',
    **regular_params)

    # **easy_params)

    # pw=500,  # IMP | (ev)
    # nbands=-50,
    # #########################################################################
    # kpts=3 * (1, ),
    # kpts=(6, 6, 6),
    # calculation='relax',
    # #########################################################################
    # parflags=pflags,
    # convergence={
    #     'energy': 1e-4,
    #     'mixing': 0.1,
    #     'maxsteps': 300,
    #     'diag': 'david',
    #     },
#__|

atoms.set_calculator(calc)

atoms.get_potential_energy()
atoms.calc.save_flev_chg("charge_den.tgz")
atoms.calc.load_flev_chg("charge_den.tgz")

# #######################################################################
import json
data_path = os.path.join(
    os.environ["repo_185b"],
    "workflow/band_path",
    "bandpath.json",
    )
with open(data_path, "r") as fle:
    data = json.load(fle)

kpts = data["kpts"]
X = data["X"]
x = data["x"]
# #######################################################################


print("")
print("atoms.calc.calc_bandstructure(kpts, atomic_projections=True)")
energies = atoms.calc.calc_bandstructure(kpts, atomic_projections=True)


if not os.path.exists("dir_bands"):
    os.makedirs("dir_bands")
with open("dir_bands/band_disp.pickle", "w") as fle:
    # pickle.dump((points, kpts, x, X, energies), fle)
    pickle.dump((kpts, x, X, energies), fle)

# COMBAK This broke when file already existed, tried a fix - 180405 - RF
shutil.move("charge_den.tgz", "dir_bands/charge_den.tgz")
atoms.write("dir_bands/out_bands.traj")



#| - __old__
#  ip = ibz_points["fcc"]
#  points = ["Gamma", "X", "W", "K", "L", "Gamma"]
#  bzpath = [ip[p] for p in points]
#  kpts, x, X = get_bandpath(bzpath, atoms.cell, npoints=300)
#__|
