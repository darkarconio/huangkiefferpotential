# NOTE: This script can be modified for different pair styles 
# See in.elastic for more info.

# Choose potential
pair_style born/coul/long 3 10.6
pair_coeff 1 1 .75946686 .303030 1.48 0 0 3
pair_coeff 1 2 .22160479 .290698 2.21 0 0 2
pair_coeff 2 2 .10335939 .256410 2.94 0 0 3
kspace_style ewald 1.0e-6

#pair_style hybrid/overlay born/coul/long 3 10.6 hk
#pair_coeff * * hk BO.hk B O
#pair_coeff 1 1 born/coul/long .75946686 .303030 1.48 0 0 3
#pair_coeff 1 2 born/coul/long .22160479 .290698 2.21 0 0 2
#pair_coeff 2 2 born/coul/long .10335939 .256410 2.94 0 0 3
#kspace_style ewald 1.0e-6

#pair_style sw
#pair_coeff * * GaN.sw Ga N

# Setup neighbor style
neighbor 1.0 nsq
neigh_modify once no every 1 delay 0 check yes 

# Setup minimization style
min_style	     cg
min_modify	     dmax ${dmax} line quadratic

# Setup output
thermo		1
thermo_style custom step temp pe press pxx pyy pzz pxy pxz pyz lx ly lz vol
thermo_modify norm no
