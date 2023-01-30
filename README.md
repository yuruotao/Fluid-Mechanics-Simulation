# Fluid-Mechanics-Simulaition

This project uses python to build a simulation of 2d fluid based on Stokes equations. This is a course project for Engineering Fluid Mechanics.

## Usage

1. Install all necessary dependencies;
2. In file named main_program.py, specify the physical parameters along with simulation parameters;
3. Run file main.py, the simulation begins

## Result

Along the process, the data would be recorded in each iteration in the form of .txt, which is in the Result folder. And after the simulation stop, the .gif file of the process would also be in the Result folder.  
For instance, use the parameters as below:

    length = 4
    width = 4
    colpts = 257
    rowpts = 257
    time = 10
    rho = 1
    mu = 0.01
    u_in = 1
    v_wall = 0
    p_out = 0

The upper parameters would yield
![The gif result](/fluid_mechanics_modeling_chn/Result/animation.gif)
And three figures  
![Fig 1](/fluid_mechanics_modeling_chn/Result/Fig1.png)
![Fig 2](/fluid_mechanics_modeling_chn/Result/Fig2.png)
![Fig 3](/fluid_mechanics_modeling_chn/Result/Fig3.png)
