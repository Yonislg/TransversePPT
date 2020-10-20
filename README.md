<h1> PPT Transeverse Code </h1>

This is a Matlab code to calculate the electron current, plasma sheath- and plasma bulk potential drop in a Pulsed Plasma Thruster channel given an input of externally appied voltage, plasma properties and PPT design parameters. This code is intended to be solved iteratively allongside an Axial code that solves the axial plasma propogation speed, density, electron temperature and magnetic field.

![PPT diagram](/DocFigs/Parallel_Plate_PPT.PNG)

<sup>source: J. K. Ziemer, Performance Scaling of Gas-Fed Pulsed Plasma Thrusters. PhD thesis, Princeton, 2001. </sup>

The PPT is the oldest flight tested electric propulsion device. It operates by igniting an arc current between two electrodes alongside a block of PTFE polymer (a.k.a. Teflon). The Teflon ablates and ionizes, forming a plasma current bridge. The magnetic field enclosed by the magnetic field propels the plasma outwards via the Lorentz force. The device operates in short pulses (~1 Hz, 5 microsreconds) which enables high instantaneous peak input power on a sattelite with limited continuous power avaliability. PPT's are considered the most simple electric propulsion devices to build. However due to the transient nature of pulsed operation, the physics are still not yet fully understood.

This code is based on a plasma model developped by Mario Merino, an associate professor and researcher at the Space Propulsion and at the PLasma and Spae Propulsion team (EP2) of Universidad Carlos III. The code itself was developped by Yonis le Grand, a graduate student of Aerospace Engineering at TU Delft and Nuclear Fusion at TU Eindhoven, while working on his graduation project under the supervision of Mario Merino.

<h2> Usage </h2>
 
To resolve the current properties of the plasma call upon the function `transversal_V3.m`, located in `+transveramodel/+main/`. This function does precalcations than calls upon `total_current.m` to resolve the model equations. It's subfunctions can be found in the folder `+transveramodel/+subfunctions/`. These include equations that calculate the collision rate, electron emissions, sheath electric field and sheath electron flux.

In the folder `+transveramodel/+testfolder/` contains functions that test the subfunctions and the whole model. In particular `test_loop.m` is a usefull script that can quickly generate potentential and current profiles for various plasma- and PPT design parameters.

<h2> Code Structure </h2>
