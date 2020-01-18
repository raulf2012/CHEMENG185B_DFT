# CHEMENG185B_DFT
Repo to store files related to running DFT simulations for CHEMENG185B

TEST commit


# Notes

## Don't have to reenter password every time you login into rice

Place the following text into this file
`$HOME/.ssh/config` (create it if needed)

```
Host rice.stanford.edu
    ControlMaster auto
    ControlPersist yes
    ControlPath ~/.ssh/%l%r@%h:%p
```
