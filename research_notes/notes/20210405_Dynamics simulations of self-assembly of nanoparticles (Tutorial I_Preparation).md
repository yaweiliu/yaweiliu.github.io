#### [HOME](../../index.html) [CONTENTS](../index.html)

## Dynamics simulations of self-assembly of nanoparticles (Tutorial I: Preparation)
 
*by Yawei Liu  @ Sydney, Australia 2021/04/04*

In following weeks, I am going to write some tutorials on dynamics simulations of the self-assembly of rod-shaped particles (or other particles with any shapes) using an overlapping-sphere model, which we did in our previous work [1,2,3]. From these tutorials, we will know:

* Using Python to create the model for a given particle (e.g. the rod-shaped particle in this case) with individual identical spheres. 
* Using Python to write a data file for the initial configurations that can be read by LAMMPS. 
* Using LAMMPS to perform dynamics simulations.
* Using Python to process the simulation data.
* Using Python to visualise the simulation system.
* AND some details you may meet in this type of simulations.

If you want to follow these tutorials, I recommend installing these softwares/packages in your computer/laptop (Linux and Mac OS are recommended) as follows:

* Python and related packages:
    * Download [Miniconda](https://docs.conda.io/en/latest/miniconda.html) for your computer/laptop and install it. This will also install Python3 in your device.
    * Install these Python packages using ```conda install package_name```: Numpy, Scipy, Pandas, Matplotlib.
    * Install Jupyterlab using ```conda install -c conda-forge jupyterlab```. (Tips: If you want to use Plato to visualise the simulation system in Jupyterlab, please install **jupyterlab2**.)

* Lammps and related packages:
    * Download [LAMMPS](https://github.com/lammps/lammps) from Github.
    * In ```/src``` folder, pre-install these packages using ```make yes-package_name```: DIPOLE, MOLECULE, RIGID, KSPACE.
    * Install LAMMPS.

* Python packages for analysis and visualisation of simulation data:
    * [Freud](https://freud.readthedocs.io/en/latest/): ```conda install -c conda-forge freud```
    * [Plato](https://plato-draw.readthedocs.io/en/latest/): ```pip install plato-draw```
    
    To use Plato in Jupyterlab, install [nodejs](https://nodejs.org/en/) and [pythreejs](https://github.com/jupyter-widgets/pythreejs) as follows:
    
    * ```conda install -c conda-forge 'nodejs>=12'```
    * ```conda install -c conda-forge pythreejs```
    * ```jupyter labextension list```
    * ```jupyter labextension install --no-build @jupyter-widgets/jupyterlab-manager```
    * ```jupyter labextension install --no-build jupyter-datawidgets/extension```
    * ```jupyter labextension install jupyter-threejs```
    * ```jupyter labextension list``````

* Also, some softwares like VMD can also be used to visualise the simulation system, I am going to write some Python scripts to convert the simulation data to the file that can be read by these softwares.


### References

[1] [Liu, Y.; Widmer-Cooper, A. A Versatile Simulation Method for Studying Phase Behavior and Dynamics in Colloidal Rod and Rod-Polymer Suspensions. J. Chem. Phys. 2019, 150 (24), 244508.](http://aip.scitation.org/doi/10.1063/1.5096193)

[2] [Sharma, A.; Wojciechowski, J. P.; Liu, Y.; Pelras, T.; Wallace, C. M.; Müllner, M.; Widmer-Cooper, A.; Thordarson, P.; Lakhwani, G. The Role of Fiber Agglomeration in Formation of Perylene-Based Fiber Networks. Cell Reports Phys. Sci. 2020, 1 (8), 100148.] (https://linkinghub.elsevier.com/retrieve/pii/S2666386420301521)

[3] [Liu, Y.; Widmer-Cooper, A. A Dissipative Particle Dynamics Model for Studying Dynamic Phenomena in Colloidal Rod Suspensions. J. Chem. Phys. 2021, 154 (10), 104120.](https://doi.org/10.1063/5.0041285)


##### Github Page / Gitee Page / Subscription
<img src="images/github_yawei.png" alt="github page" width="80" height="80" />
<img src="images/gitee_yawei.png" alt="gitee page" width="80" height="80" />
<img src="images/wechat.png" alt="wechat" width="80" height="80" />

<p>&copy; 2021 Yawei Liu. All content licensed under the <a href="https://creativecommons.org/licenses/by-nc/4.0/legalcode#languages">Creative Commons Attribution-NonCommercial License 4.0 International (CC BY-NC 4.0)</a>.</p>

--
#### [HOME](../../index.html) [CONTENTS](../index.html)
