#### [HOME](../../index.html) [CONTENTS](../index.html)

## Dynamics simulations of self-assembly of nanoparticles (Tutorial V: Creating LAMMPS script)
 
*by Yawei Liu  @ Sydney, Australia 2021/05/02*

In this tutorial, I am going to show a LAMMPS script to run a compression Langevin simulation for a group of rod-shaped nanoparticles.

In our model, rod spheres interact with other spheres via a continuous pseudo-hard-core potential, in the form of $U=50(50/49)^{49} \epsilon [(\sigma/r)^{50} -(\sigma/r)^{49}]$ truncated and shifted at $r^{\alpha\beta}_{cut}=(50/49)\sigma$. Here, $r$ is the center-to-center distance between the spheres, $\epsilon$ is the energy parameter, $\sigma$ is the distance parameter [1,2]. During the Langevin simulations, all spheres are subjected to three forces: a conservative force $F^C$ computed via the pairwise interactions described above; a friction force $F^F=-(m/\gamma)v$ with $m$ the sphere mass, $\gamma$ the damping factor, and $v$ the sphere velocity; and a random force $F^R\propto\sqrt{k_BT m/(\Delta t\gamma)}$ with $k_B$ the Boltzmann constant, $T$ the desired temperature, and $\Delta t$ the time step. Here, we set the mass of each rod sphere to be the same, i.e. $m_p=m$. All simulations are carried out at a dimensionless temperature $k_BT/\epsilon=1$. The velocity-Verlet algorithm was used to integrate the equations of motion with a time step $\Delta t=0.005\tau$ where $\tau=D\sqrt{m/(k_BT)}$. A Berendsen barostat with a time constant of $5\tau$ was applied to control the pressure.

Here is the LAMMPS script:

```
variable    rodL     equal 5.0
variable    mypress  equal 1.0
variable    myPstart equal 1.0
variable    myPend   equal 3.0
variable    mytemp   equal 1.0
variable    myrand   equal 4003

#####0 relax  1 compress
variable    restart_flag equal 0 

units       lj
dimension   3
boundary    p p p
atom_style  hybrid bond dipole

if "${restart_flag} > 0" then &
    "read_data restart.dat" &
    "jump input_eq.lmp main"

read_data  init.dat
mass       * 1
velocity   all create ${mytemp} ${myrand} rot yes dist gaussian
##########################################################
label      main

###potential
pair_style mie/cut 1.02040816327
pair_modify shift yes
pair_coeff   *   *  1.0 1.0 50 49 

timestep 0.005

group origin type 1
group solid type 1 2 3
group output type 1

variable commcutoff equal ${rodL}/2.+0.1
comm_modify mode single cutoff ${commcutoff}
neigh_modify exclude molecule/intra solid

fix rigid solid rigid/small molecule
if "${restart_flag} ==0 " then &
   "fix press all press/berendsen iso ${mypress} ${mypress} 5.0"
if "${restart_flag} >0 " then &
   "fix press all press/berendsen aniso ${myPstart} ${myPend} 5.0"
fix langevin all langevin 1.0 1.0 1 ${myrand} zero yes

variable step equal step
variable time equal time
variable lx equal lx
variable ly equal ly
variable lz equal lz
variable Nc equal count(origin)

variable temp equal temp
variable pxx equal pxx
variable pyy equal pyy
variable pzz equal pzz

thermo_style custom step atoms temp press pxx lx ly lz v_Nc
thermo_modify flush yes
thermo 1000

fix log all print 1000 "${step} ${time} ${temp} ${pxx} ${pyy} ${pzz} ${lx} ${ly} ${lz} ${Nc}" title "step time temp pxx pyy pzz lx ly lz Nr" file result_thermo.log screen no

dump atom output custom 5000 result_atoms.log id mol type x y z mux muy muz
dump_modify atom sort id

run 1000000
write_data  restart.dat.* nofix nocoeff

```

When ```restart_flag==0```, LAMMPS will read the initial configuration file ```init.dat``` and equilibrate the system at a low pressure (thus an isotropic phase will be otbained). The final configuration will be saved as ```restart.dat```.

When ```restart_flag==1```, LAMMPS will read the file ```restart.dat``` and compress the system from a low pressure to a high pressure, and we will observe the isotropic-nematic-smectic-crystal phase transitions as the pressure increases.

### References

[1] [Liu, Y.; Widmer-Cooper, A. A Versatile Simulation Method for Studying Phase Behavior and Dynamics in Colloidal Rod and Rod-Polymer Suspensions. J. Chem. Phys. 2019, 150 (24), 244508.](http://aip.scitation.org/doi/10.1063/1.5096193)

[2] [Liu, Y.; Widmer-Cooper, A. A Dissipative Particle Dynamics Model for Studying Dynamic Phenomena in Colloidal Rod Suspensions. J. Chem. Phys. 2021, 154 (10), 104120.](https://aip.scitation.org/doi/10.1063/5.0041285)


##### Github Page / Gitee Page / Subscription
<img src="images/github_yawei.png" alt="github page" width="80" height="80" />
<img src="images/gitee_yawei.png" alt="gitee page" width="80" height="80" />
<img src="images/wechat.png" alt="wechat" width="80" height="80" />

<p>&copy; 2021 Yawei Liu. All content licensed under the <a href="https://creativecommons.org/licenses/by-nc/4.0/legalcode#languages">Creative Commons Attribution-NonCommercial License 4.0 International (CC BY-NC 4.0)</a>.</p>

--
#### [HOME](../../index.html) [CONTENTS](../index.html)
