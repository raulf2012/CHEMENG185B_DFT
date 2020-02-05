#!/usr/bin/env python
#
#################
#SBATCH --job-name=test00
#################
#SBATCH --time=15:00
#################
#SBATCH --ntasks=1
#################
#SBATCH --cpus-per-task=1
#################
#SBATCH --nodes=1
#################
#SBATCH --ntasks-per-node=16
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

# "nodes": "1",  # --nodes


"""Simple DFT optimization example.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
from ase import io

from espresso.espresso import Espresso
from espresso.espresso import SiteConfig
#__|


# Read atoms object
try:
    atoms = io.read("init.traj")
except:
    atoms = io.read("init.cif")

#| - Parallization
#  ni = 1
#  nk = 6
#  nt = 2
#  nd = 16

ni = 1
nk = 4
nt = 1
nd = 4

pflags = '-ni ' + str(ni) + ' -nk ' + str(nk) + ' -nt ' + str(nt)+ ' -nd ' + str(nd)
#__|


#| - QE Calculator
site_config = SiteConfig(
    # scheduler=None,
    # usehostfile=False,
    scratchenv='SCRATCH',
    )

calc_new = Espresso(
    atoms=None,
    pw=350.0,
    dw=None,
    fw=None,
    nbands=-10,
    # #########################################################################
    kpts=3 * (3, ),

    kptshift=(0, 0, 0),
    fft_grid=None,
    calculation='relax',
    ion_dynamics='bfgs',
    nstep=None,
    constr_tol=None,
    fmax=0.05,
    cell_dynamics=None,
    press=None,
    dpress=None,
    cell_factor=None,
    cell_dofree=None,
    dontcalcforces=False,
    nosym=False,
    noinv=False,
    nosym_evc=False,
    no_t_rev=False,
    xc='PBE',
    beefensemble=False,
    printensemble=False,
    psppath=None,
    spinpol=False,
    noncollinear=False,
    spinorbit=False,
    # #########################################################################
    outdir="calcdir",

    txt=None,
    calcstress=False,
    smearing='fd',
    sigma=0.1,
    fix_magmom=False,
    isolated=None,
    U=None,
    J=None,
    U_alpha=None,
    U_projection_type='atomic',
    nqx1=None,
    nqx2=None,
    nqx3=None,
    exx_fraction=None,
    screening_parameter=None,
    exxdiv_treatment=None,
    ecutvcut=None,
    tot_charge=None,
    charge=None,
    tot_magnetization=-1,
    occupations='smearing',
    dipole={'status': False},
    field={'status': False},
    output={
        # Default
        # 'disk_io': 'default',
        # 'avoidio': False,
        # 'removewf': True,
        # 'removesave': False,
        # 'wf_collect': False,

        'disk_io': 'default',
        'avoidio': True,
        'removewf': True,
        'removesave': True,
        'wf_collect': False,
        },
    convergence={
        # Default
        # 'energy': 1e-6,
        # 'mixing': 0.7,
        # 'maxsteps': 100,
        # 'diag': 'david',

        'energy': 1e-3,
        'mixing': 0.1,
        'maxsteps': 300,
        'diag': 'david',
        },
    startingpot=None,
    startingwfc=None,
    ion_positions=None,
    # #########################################################################
    #  parflags=None,
    parflags=pflags,

    alwayscreatenewarrayforforces=True,
    verbose='low',
    # automatically generated list of parameters
    # some coincide with ase-style names
    iprint=None,
    tstress=None,
    tprnfor=None,
    dt=None,
    lkpoint_dir=None,
    max_seconds=None,
    etot_conv_thr=None,
    forc_conv_thr=None,
    tefield=None,
    dipfield=None,
    lelfield=None,
    nberrycyc=None,
    lorbm=None,
    lberry=None,
    gdir=None,
    nppstr=None,
    nbnd=None,
    ecutwfc=None,
    ecutrho=None,
    ecutfock=None,
    force_symmorphic=None,
    use_all_frac=None,
    one_atom_occupations=None,
    starting_spin_angle=None,
    degauss=None,
    nspin=None,
    ecfixed=None,
    qcutz=None,
    q2sigma=None,
    x_gamma_extrapolation=None,
    lda_plus_u=None,
    lda_plus_u_kind=None,
    edir=None,
    emaxpos=None,
    eopreg=None,
    eamp=None,
    clambda=None,
    report=None,
    lspinorb=None,
    esm_w=None,
    esm_efield=None,
    esm_nfit=None,
    london=None,
    london_s6=None,
    london_rcut=None,
    xdm=None,
    xdm_a1=None,
    xdm_a2=None,
    electron_maxstep=None,
    scf_must_converge=None,
    conv_thr=None,
    adaptive_thr=None,
    conv_thr_init=None,
    conv_thr_multi=None,
    mixing_beta=None,
    mixing_ndim=None,
    mixing_fixed_ns=None,
    ortho_para=None,
    diago_thr_init=None,
    diago_cg_maxiter=None,
    diago_david_ndim=None,
    diago_full_acc=None,
    efield=None,
    tqr=None,
    remove_rigid_rot=None,
    tempw=None,
    tolp=None,
    delta_t=None,
    nraise=None,
    refold_pos=None,
    upscale=None,
    bfgs_ndim=None,
    trust_radius_max=None,
    trust_radius_min=None,
    trust_radius_ini=None,
    w_1=None,
    w_2=None,
    wmass=None,
    press_conv_thr=None,

    # site=None,
    site=site_config,
    )
# __|

atoms.set_calculator(calc_new)

atoms.get_potential_energy()

# atoms.write("final.traj")

# Read atoms object from qe log file
# io.read("calcdir/log")

