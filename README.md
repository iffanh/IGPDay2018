# IGPDay2018
An LBM code provided for people to see, test, and play around during IGP Day 29/10/2018. 

This code was made for the TPG4560 - Petroleum Engineering, Specialization Project as a means to understand the basic concept of Lattice Boltzmann Modeling.

In this repository the codes for poiseuille flow case and flow around a cylinder case are presented on Python platform. These two cases are separated by the folders, named after the cases themselves. To start the simulation, run the main file "lbm_cylinder.py" and "lbm_poiseuille.py". The main input data can be changed in the corresponding input files "lbm_inputdata_cylinder.py" and "lbm_poiseuille.py".

You will probably need to install seaborn, a very nice package to generate heatmap easily. To install, type

pip install seaborn

The code will run for 5 to 6 hours, depends on the computer capabilities. The result will be saved as .png image files for every 100 time steps. 

Bear in mind that the code is far from perfect and thus when the input does not match certain requirements, the result may be unphysical due to LBM limitation for compressible flow. 

Any critics or suggestion can be sent to my email, muhammih@stud.ntnu.no or iffan.hannanu@yahoo.com or iffan.hannanu@gmail.com. I will be very happy to see some suggestions!

by Muhammad Iffan Hannanu
Saturday, 20/10/2018, 18:59
In Schlumberger Computer Lab, 4th floor
Department of Geoscience and Petroleum
Lerkendal
Trondheim, Norway