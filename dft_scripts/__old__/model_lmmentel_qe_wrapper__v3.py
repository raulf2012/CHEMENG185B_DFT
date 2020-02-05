#!/usr/bin/env python
#| - SLURM HEADER
#################
#SBATCH --job-name=test00
#################
#SBATCH --time=60:00
#################
#SBATCH --ntasks=8
#################
#SBATCH --nodes=1
#################
#SBATCH --mem-per-cpu=2G
#################
#get emailed about job BEGIN, END, and FAIL
#SBATCH --mail-type=END,FAIL
#################
#who to send email to; please change to your email
#SBATCH  --mail-user=flores12@stanford.edu
#################
#a file for job output, you can check job progress
#SBATCH --output=job.out
#################
# a file for errors from the job
#SBATCH --error=job.err
#################

#################
# #SBATCH --cpus-per-task=1
#################
# #SBATCH --ntasks-per-node=1
#__|


"""Simple DFT optimization example.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

from ase import io

from espresso import Espresso
from espresso import SiteConfig
#__|


#| - Read atoms object
try:
    atoms = io.read("init.traj")
except:
    atoms = io.read("init.cif")
#__|

#| - Parallization
#  ni = 1
#  nk = 6
#  nt = 2
#  nd = 16

#  ni = 1
#  nk = 4
#  nt = 1
#  nd = 4

ni = 1
nk = 1
nt = 1
nd = 1

pflags = '-ni ' + str(ni) + ' -nk ' + str(nk) + ' -nt ' + str(nt)+ ' -nd ' + str(nd)
#__|

#| - QE Calculator

site_config = SiteConfig(
    scheduler="SLURM",
    # usehostfile=False,
    scratchenv='SCRATCH',
    )

calc = Espresso(
    pw=500,  # IMP | (ev)
    nbands=-50,
    # #########################################################################
    # kpts=3 * (1, ),
    kpts=(6, 6, 6),
    kptshift=(0, 0, 0),
    calculation='relax',
    ion_dynamics='bfgs',
    fmax=0.05,
    xc='PBE',  # IMP
    # #########################################################################
    outdir="calcdir",
    U_projection_type='atomic',
    site=site_config,
    tot_magnetization=-1,
    occupations='smearing',
    dipole={'status': False},
    field={'status': False},
    output={
        'disk_io': 'default',
        'avoidio': True,
        'removewf': True,
        'removesave': True,
        'wf_collect': False,
        },
    convergence={
        'energy': 1e-4,
        'mixing': 0.1,
        'maxsteps': 300,
        'diag': 'david',
        },
    smearing='fd',
    sigma=0.1,
    # #########################################################################
    # parflags=pflags,
    alwayscreatenewarrayforforces=True,
    verbose='low',
    )

# __|

atoms.set_calculator(calc)

try:
    print("calc.localtmp:", calc.localtmp)
except:
    print("sidjfisjifjisdji")

atoms.get_potential_energy()



#| - __old__
#  print(40 * "*")
#  print(atoms.calc.site)
# Read atoms object from qe log file
# io.read("calcdir/log")
#  vossjo_espresso_params = dict(
#      #  mode="bfgs",
#      mode="relax",
#      opt_algorithm="bfgs",
#      )
#  calc = Espresso(
#      pw=350.0,
#      kpts=3 * (1, ),
#      **vossjo_espresso_params)
# #####################################
# atoms.write("final.traj")

#  print("SLURM_JOB_ID:", os.environ["SLURM_JOB_ID"])
#  print("IDJIFJSDF")
#  print(os.environ)
# print("SUBMITDIR:", os.environ["SUBMITDIR"])
#  print("SUBMITDIR:", os.environ["SUBMITDIR"])
#  SLURM_SUBMIT_DIR
#__|
