#### [HOME](../../index.html) [CONTENTS](../index.html)

## Units conversion in dissipative particle dynamics (DPD)

*by Yawei Liu  @Sydney, Australia 2020/01/16*

### 1. DPD model

In the DPD model [1], one DPD bead represents $N_m$ fluid (e.g. water) molecules. Here, $N_m$ is also called coarse-graining (CG) degree. Each pair of DPD beads $i$ and $j$ interact with each other via a conservative force $\textbf{F}_{ij}^C = A (1-r_{ij}/r_c) \textbf{e}_{ij}$, a dissipative force
$\textbf{F}_{ij}^D = -\gamma(1-r_{ij}/r_c)^2 (\textbf{e}_{ij} \cdot \textbf{v}_{ij}) \textbf{e}_{ij}$, and a random force $\textbf{F}_{ij}^R = \xi (1-r_{ij}/r_c) \theta_{ij} (\Delta t)^{-1/2} \textbf{e}_{ij}$ and all forces vanish at $r\ge r_c$. Here,

* $r_{ij}=|\textbf{r}_{ij}|$: the centre-to-centre distance between the two beads.
* $\textbf{e}_{ij}=\textbf{r}_{ij}/r_{ij}$: the unit vector pointing between the two beads.
* $\textbf{v}_{ij}$: the vector difference in velocities between the two beads.
* $A$: the repulsion coefficient.
* $\gamma$: the dissipative coefficient.
* $\xi$: the the noise strength, and $\xi=\sqrt{2k_BT\gamma}$ with $k_B$ the Boltzmann constant and $T$ the temperature.
* $\theta_{ij}$: a Gaussian white noise variable.
* $\Delta t$: the simulation time step.
* $r_c$: the cutoff distance, and also can be treated as the size (diameter) of the DPD bead.

### 2. Determining paramters

* $r_c$, $k_BT$, $m$ (bead mass) and $\tau$ (time)

In simulations, the reduced units are often used. For DPD model, all units are often scaled by the length unit $r_c$, the mass unit $m$, the energy unit $k_BT$, and the time unit $\tau$. Hence, $r_c^* = m^* = (k_BT)^* =\tau^* =1$ in the simulations. The superscript asterisk ($^*$) means the quality is in reduced units. Units of other quantities are from these basic units. For example, the unit for the mass density is $m/r_c^3$, the unit for the diffusion constant is $r_c^2/\tau$.

* $\rho$ (density)

The DPD simulations are normally carried out within a $NVT$ ensemble, and the number of particles are often determined by setting $\rho^*=3$. When $\rho^*>2$, a simple scaling relation between the density and excess pressure exits [1]. In principle the density chosen for the simulation is a free parameter, but for efficiency reasons one would thus choose the lowest possible density where the scaling relation still holds.

* $A$ and $N_m$

In order to match the compressibility of DPD fluid with a liquid having the dimensionless compressibility of $\Psi$, the interparticle repulsion coefficient $A$ is often chosen through
$$A\approx \frac{\Psi N_m -1}{0.2\rho^*}(k_BT)^*.$$
For water, $\Psi\approx 16$ yields $A^*=25$ for $N_m=1$ and $A^*=104$ for $N_m=4$.

* $\gamma$ and $\Delta t$

As a reasonable compromise between fast temperature equilibration, a fast simulation and a stable, physically meaningful system, simulation with $\Delta t^*=0.04$ and $\gamma^*=4.5$ is often recommended.
However, above choice for $\Delta t$ and $\gamma$ yields very small Schmidt number (i.e. the ratio of viscosity to diffusion) ($Sc\sim1$) compared to the real fluid such as water ($Sc\sim10^3$). A possible solution to this problem is increasing $\gamma$. At the same time though, $\Delta t$ would have to be reduced to maintain the temperature control.

* **There is no unique way to determine the parameters for DPD model.** For a given physical problem, with a characteristic length scale, we may always put a given number of DPD particles and parametrize the model in order to recover some macroscopic information (e.g. compressibilities, viscosity and diffusion constant).

### 3. Units conversion
The conversation between DPD units ($r_c$, $k_BT$, $m$ and $\tau$) and real units ([m], [J], [kg] and [s]) is obtained from the key macroscopic information recovered by the model. If quantities in DPD units are labeled with $_{sim}$ and in real units are labeled with $_{real}$, one would have

* $m=N_m M_{real}/N_A$ with $M_{real}$ the mole mass of the fluid (e.g. water) and $N_A$ the Avogadro constant.
* $r_c=(m\rho_{sim}/\rho_{real})^{1/3}$ with $\rho_{real}$ the mass density of the fluid.
* $k_BT = k_B \cdot T$ with $k_B=1.38064853e-23$ J$\cdot$K and $T$ is the system temperature in K.
* The time unit $\tau$ can be chosen in different ways.
  * By taking the long-term self-diffusion constant into account, $\tau=N_m r_c^2 D_{sim}/D_{real}$ with $D$ the diffusion constant. Sometimes, this relation is also used to determine $N_m$.
  * If the viscous processes are the main parts, then $\tau=r_c^2 \nu_{sim}/\nu_{real}$ with $\nu$ the kinematic viscosity.
  * Or $\tau=r_c/U_{ref}$ in which the reference velocity $U_{ref}$ is chosen as either system thermal velocity $(k_BT/m)^{1/2}$ or the characteristic velocity of the real flow.

### 4. Example
An $NVT$ simulation with $1000$ DPD beads in a cubic box is carried out by LAMMPS. The parameters are: $\rho^*=3$; $A=104$ for $N_m=4$; $\gamma=100$; $\Delta t=0.005$. Then, comparing the DPD fluid with the water at $T=298$ K, we have:

* Mass $m=1.1969\times10^{-25}$ kg ($m=N_m M_{real}/N_A$)
* Length $r_c=7.1125\times10^{-10}$ m ($r_c=(m\rho_{sim}/\rho_{real})^{1/3}$)
* Energy $k_BT = 4.1143\times10^{-21}$ J ($k_BT = k_B \cdot T$)
* Time $\tau = 3.8363\times10^{-12}$ s ($\tau=r_c\sqrt{m/k_BT}$ )

As a results, the kinematic viscosity for the DPD fluid is $\nu=1.569 r_c^2/\tau = 2.069\times10^{-7}$ m$^2$/s ($\nu=8.935\times10^{-7}$ m$^2$/s for water).

Also see our published papers for some particular examples in charged systems [2,3] and vapour-liquid systems [4].

### References
[1] [Groot, R. D.; Warren, P. B. Dissipative Particle Dynamics: Bridging the Gap between Atomistic and Mesoscopic Simulation. J. Chem. Phys. 1997, 107 (11), 4423–4435.](http://aip.scitation.org/doi/10.1063/1.474784)

[2] [Zhang, H.; Liu, Y.; Shahidan, M. F. S.; Kinnear, C.; Maasoumi, F.; Cadusch, J.; Akinoglu, E. M.; James, T. D.; Widmer‐Cooper, A.; Roberts, A.; et al. Direct Assembly of Vertically Oriented, Gold Nanorod Arrays. Adv. Funct. Mater. 2021, 31 (6), 2006753.](https://onlinelibrary.wiley.com/doi/10.1002/adfm.202006753)

[3] [Wei, J.; Liu, Y.; Song, F. Coarse-Grained Simulation of the Translational and Rotational Diffusion of Globular Proteins by Dissipative Particle Dynamics. J. Chem. Phys. 2020, 153 (23), 234902.](http://aip.scitation.org/doi/10.1063/5.0025620)

[4] [Liu, Y.; Bernardi, S.; Widmer-Cooper, A. Stability of Pinned Surface Nanobubbles against Expansion: Insights from Theory and Simulation. J. Chem. Phys. 2020, 153 (2), 024704.](http://aip.scitation.org/doi/10.1063/5.0013223)

##### Github Page / Gitee Page / Subscription
<img src="images/github_yawei.png" alt="github page" width="80" height="80" />
<img src="images/gitee_yawei.png" alt="gitee page" width="80" height="80" />
<img src="images/wechat.png" alt="wechat" width="80" height="80" />

<footer>
    <script async src="//busuanzi.ibruce.info/busuanzi/2.3/busuanzi.pure.mini.js"></script>
    <span id="busuanzi_container_page_pv" style='display:none'>
      <h6>view <span id="busuanzi_value_page_pv">       </span> times</h6>
    </span>
</footer>

<p>&copy; 2021 Yawei Liu. All rights reserved.</p>

--
#### [HOME](../../index.html) [CONTENTS](../index.html)

