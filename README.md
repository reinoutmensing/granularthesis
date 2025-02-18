# granularthesis

To execute the simulation, you need to run three codes in this order: 
 - CSCWalls.cpp
 - CSCInit.cpp
 - CSCRun.cpp

## CSCWalls
This code is used to define the particle and material properties and the geometry:

We use particles that are non-dimensionalised such that their diameter d=1 and mass m=1. Gravitational acceleration is also non-dimensionalised to g=(0,0,1); thus, the gravitational time-scale is t_g=sqrt(d/g)=1. We use a linear spring-damper contact law with collision time t_c=t_g/20, restitution r=0.88, and a sliding friction of mu=0.5. 

This code defines a box of size W=30 in x-, L=20 in y-, and H=30 in z-direction. The walls consist of flat walls with particles attached to the wall to make it more frictional.

## CSCInit

This code restarts from the box created in CSCWalls, and fills the box with particles

## CSCRun

This code restarts from the box created in CSCInit, and adds a shear velocity by moving the walls at a steady velocity in y-direction: The right wall and the right half of the bottom wall is moved at a velocity (0,v,0), the left wall and the left half of the bottom wall is moved at a velocity (0,-v,0).

## readMercury.CG 

This script reads and processes coarse-grained simulation data from MercuryDPM .stat files. It extracts metadata from the file headers, determines the spatial dimensions, and organizes variables into structured fields. The script reconstructs coordinates (x, y, z), physical quantities (such as density, momentum, and stress tensors), and statistical properties of particle sizes. Additionally, it derives secondary variables, such as velocity fields and kinetic stress, by computing relations from the available data. The output is a structured data object containing all relevant fields, making it easier to analyze and visualize the simulation results.

## Author

These four files have been written by Thomas Weinhart. However, they have been adjusted.

## steadystate.m

After running ./MercuryCG CSCRun, this file is used to see when the steadystate kicks in. After finding this steady state ./MercuryCG CSCRun -coordinates XZ -n 100 -tmin 400 -timeaverage -w 3 (this is just an example) should be reran with -tmin anywhere after the steady state has been reached.

## main.m

This script processes simulation data from MercuryDPM .stat files and analyzes shear band behavior in granular materials. It reads and preprocesses the data, computes shear stress and shear rate, and performs curve fitting to characterize the velocity profile and inertial number variations. The script generates multiple contour plots, shear rate vs. shear stress comparisons, and statistical analyses of key parameters. Additionally, it explores the relationship between shear stress, shear rate, and inertial number through various fitting models and plots. The results include visualizations of shear band width, height, and center.

## Janssen.m

It extracts contact stress components (\(\sigma_{xx}, \sigma_{yy}, \sigma_{zz}\)), applies NaN handling to ensure data consistency, and generates contour plots to visualize stress variations. Additionally, it computes gravitational stress (\(\sigma_L\)), evaluates wall traction (\(T_z\)), and examines the relationship between stress and height. The script also performs cumulative integration to derive gravitational stress and compares it with the measured stress field.

## shearpath.m

shearpath.m reads and sorts all available `.stat` files, extracts unique height values, and identifies the top-filled height in the system. The script then scans through the vertical axis to detect zero crossings in velocity data and computes key properties such as the average x-position at high z-levels and the first height where this value is exceeded. The results include plots visualizing shear band behavior, with a focus on detecting variations across different heights. Additionally, it provides insights into shear band fluctuations by tracking height-dependent changes in velocity field zero crossings.

## critical.m

The script calculates properties such as the average x-position at high z-levels and the first height where this value is exceeded. Additionally, it classifies files based on predefined categories and assigns colors accordingly. 

## contour.m

The script detects zero crossings in the velocity field along the vertical axis and calculates the average x-position at high z-levels, as well as the first height where this value is exceeded. Additionally, it evaluates contour values, accounts for oscillatory contributions, and normalizes them based on system width. The final outputs include visualizations of shear band fluctuations, cumulative contour values, and oscillation contributions, providing insights into granular flow dynamics.

## stat_specific

stat specific contains the files that are used for the experiment. However, since github does not allow big files, the files contained here are only a couple of smaller files to show the idea.

## figures

Figures contains the figures. These are now the figures that have been created with these smaller files that can be uploaded to github.
