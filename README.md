# CHEMENG185B_DFT
Repo to store files related to running DFT simulations for CHEMENG185B

TEST commit


# Notes


## Getting QE up and running on rice

### Getting OpenMPI loaded properly

source /etc/profile.d/lmod.sh
module load openmpi/3.0.0

### LD_LIBRARY_PATH on Sherlock (working fine with new qe wrapper)
Printing following environmental variable: LD_LIBRARY_PATH
-----------------------------------------------------
/home/users/vossj/suncat/s2/ompi2.1.0/lib
/share/software/user/restricted/icc/2017.u2/lib/intel64
/share/software/user/restricted/ifort/2017.u2/lib/intel64
/share/software/user/restricted/imkl/2017.u2/lib/intel64
/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/intel64
/opt/intel/composer_xe_2013_sp1.1.106/mkl/lib/intel64
/home/users/vossj/suncat/lib
/home/users/vossj/suncat/lib64
/home/users/vossj/suncat/lib/s2.0
/share/software/user/restricted/imkl/2017.u2/mkl/lib/intel64
/share/software/user/restricted/ifort/2017.u2/lib/intel64
/share/software/user/restricted/icc/2017.u2/lib/intel64




## Don't have to reenter password every time you login into rice

Place the following text into this file
`$HOME/.ssh/config` (create it if needed)

```
Host rice.stanford.edu
    ControlMaster auto
    ControlPersist yes
    ControlPath ~/.ssh/%l%r@%h:%p
```

# Bi pseudopotential file
https://www.quantum-espresso.org/upf_files/Bi.pbe-dn-rrkjus_psl.1.0.0.UPF

I also put it in the Team Drive folder
