# FSL fork of MP-SPDZ

To sync the fork, run the following:
```zsh
git fetch upstream && git merge --no-commit --no-ff upstream/master
```
Fix any merge conflicts, and then commit+push your changes. 

This assumes you have already added the original MP-SPDZ repository as an
upstream remote:

```zsh
git remote add upstream https://github.com/data61/MP-SPDZ.git
```