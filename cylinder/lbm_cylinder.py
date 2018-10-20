#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
from numpy import *
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import lbm_inputdata_cylinder as ini             #import from theinput file
import os

####################################################### SIMULATION ###########################################################################

for t in range(ini.T):

    #For smoothing, test friday 19/10/2018
    if t == 0:
        ini.tau0 = 2.*ini.tau0

    if t == 100:
        ini.tau0 = ini.tau0/2.

    #Zou and He velocity BCs on west side
    ini.rho[1,:] = (ini.f[1,:,0] + ini.f[1,:,2] + ini.f[1,:,4] + 2.*(ini.f[1,:,3] + ini.f[1,:,7] + ini.f[1,:,6])) / (1 - ini.ux0*(1 + 1e-4*sin(2.*arange(ini.sizeY_+2)*pi/ini.sizeY_)))
    ru = ini.rho[1,:]*ini.ux0*(1 + 1e-4*sin(2*arange(ini.sizeY_+2)*pi/ini.sizeY_ + 2))
    ini.f[1,:,1] = abs(ini.f[1,:,3] + (2./3.)*ru)
    ini.f[1,:,5] = abs(ini.f[1,:,7] + (1./6.)*ru - (1./2.)*(ini.f[1,:,2] - ini.f[1,:,4]))# + sin(t*pi/7)/2
    ini.f[1,:,8] = abs(ini.f[1,:,6] + (1./6.)*ru - (1./2.)*(ini.f[1,:,4] - ini.f[1,:,2]))# - sin(t*pi/7)/2

    # ... computing density for imaging
    ini.rho[:,:] = 0.
    for a in range(9):
        ini.f[:,:,a] = np.where(ini.f[:,:,a] > ini.f_tol, ini.f[:,:,a], ini.f_tol)
        ini.rho[:,:] += np.where(ini.m[:,:] == 0, ini.f[:,:,a], 0)

    #Streaming step
    ini.ftemp[:,:,0] = ini.f[:,:,0]
    for j in range(ini.sizeY_+2):
        j_n = (j-1) if j > 0 else (ini.sizeY_+1)
        j_p = (j+1) if j < (ini.sizeY_ + 1) else 0 
        for i in range(1,ini.sizeX_+ 1):
            i_n = (i-1) 
            i_p = (i+1)
            if (ini.m[i,j] <> 1):
                #For streaming part, if adjacent grid is a m, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction        
                if (ini.m[i_p,j] <> 1):    ini.ftemp[i_p,j,1] = ini.f[i,j,1]
                else:                      ini.ftemp[i,j,3] = ini.f[i,j,1]

                if (ini.m[i,j_p] <> 1):    ini.ftemp[i][j_p][2] = ini.f[i][j][2]
                else:                      ini.ftemp[i][j][4] = ini.f[i][j][2]

                if (ini.m[i_n,j] <> 1):    ini.ftemp[i_n,j,3] = ini.f[i,j,3]
                else:                      ini.ftemp[i,j,1] = ini.f[i,j,3]

                if (ini.m[i,j_n] <> 1):    ini.ftemp[i,j_n,4] = ini.f[i,j,4]
                else:                      ini.ftemp[i,j,2] = ini.f[i,j,4]

                if (ini.m[i_p,j_p] <> 1):  ini.ftemp[i_p,j_p,5] = ini.f[i,j,5]
                else:                      ini.ftemp[i,j,7] = ini.f[i,j,5]

                if (ini.m[i_n,j_p] <> 1):  ini.ftemp[i_n,j_p,6] = ini.f[i,j,6]
                else:                      ini.ftemp[i,j,8] = ini.f[i,j,6]

                if (ini.m[i_n,j_n] <> 1):  ini.ftemp[i_n,j_n,7] = ini.f[i,j,7]
                else:                      ini.ftemp[i,j,5] = ini.f[i,j,7]

                if (ini.m[i_p,j_n] <> 1):  ini.ftemp[i_p,j_n,8] = ini.f[i,j,8]
                else:                      ini.ftemp[i,j,6] = ini.f[i,j,8]

    # ... and then computing macroscopic density and velocity for each lattice point, after shifting
    ini.rho[:,:] = 0.
    ini.ux[:,:] = 0.
    ini.uy[:,:] = 0.
    for a in range(9):
        ini.rho[1:ini.sizeX_+1,:] += ini.ftemp[1:ini.sizeX_+1,:,a]    
        ini.ux[1:ini.sizeX_+1,:] += ini.e_[0,a]*ini.ftemp[1:ini.sizeX_+1,:,a]
        ini.uy[1:ini.sizeX_+1,:] += ini.e_[1,a]*ini.ftemp[1:ini.sizeX_+1,:,a]
    ini.ux[1:ini.sizeX_+1,:] = ini.ux[1:ini.sizeX_+1,:]/ini.rho[1:ini.sizeX_+1,:]
    ini.uy[1:ini.sizeX_+1,:] = ini.uy[1:ini.sizeX_+1,:]/ini.rho[1:ini.sizeX_+1,:]
    ini.u[1:ini.sizeX_+1,:] = sqrt((ini.ux[1:ini.sizeX_+1,:]**2 + ini.uy[1:ini.sizeX_+1,:]**2)/2)

    #Calculating drag force on the cylinder
    # F_drag = 0
    # for j in range(2,ini.sizeY_):
    #     j_n = (j-1) 
    #     j_p = (j+1)
    #     for i in range(2,ini.sizeX_):
    #         i_n = (i-1) 
    #         i_p = (i+1)
    #         if (ini.m[i,j] == 0): 
    #             if (ini.m[i_p,j] == 1):    F_drag += ini.ux[i,j]**2
    #             if (ini.m[i,j_p] == 1):    F_drag += ini.uy[i,j]**2
    #             if (ini.m[i_n,j] == 1):    F_drag += ini.ux[i,j]**2
    #             if (ini.m[i,j_n] == 1):    F_drag += ini.uy[i,j]**2
    #             if (ini.m[i_p,j_p] == 1):  F_drag += ini.u[i,j]**2
    #             if (ini.m[i_n,j_p] == 1):  F_drag += ini.u[i,j]**2
    #             if (ini.m[i_n,j_n] == 1):  F_drag += ini.u[i,j]**2
    #             if (ini.m[i_p,j_n] == 1):  F_drag += ini.u[i,j]**2
    # #F_drag = F_drag
    # ini.Ft[t] = F_drag 

    fct1 = ini.w[0]*ini.rho[:,:]
    fct2 = ini.w[1]*ini.rho[:,:]
    fct3 = ini.w[2]*ini.rho[:,:]

    ini.uxeq[:,:] = ini.ux[:,:]                                       #uxeq will incorporate external forces, if any
    ini.uyeq[:,:] = ini.uy[:,:] 

    ini.uxsq[:,:] = ini.uxeq[:,:]*ini.uxeq[:,:]
    ini.uysq[:,:] = ini.uyeq[:,:]*ini.uyeq[:,:]

    ini.uxuy5[:,:] = ini.uxeq[:,:] + ini.uyeq[:,:]
    ini.uxuy6[:,:] = -ini.uxeq[:,:] + ini.uyeq[:,:]
    ini.uxuy7[:,:] = -ini.uxeq[:,:] - ini.uyeq[:,:]
    ini.uxuy8[:,:] = ini.uxeq[:,:] - ini.uyeq[:,:]

    ini.usq[:,:] = ini.uxsq[:,:] + ini.uysq[:,:]

    ini.feq[:,:,0] = fct1*(1.                                                   - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,1] = fct2*(1. + ini.c_eq[0]*ini.uxeq[:,:] + ini.c_eq[1]*ini.uxsq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,2] = fct2*(1. + ini.c_eq[0]*ini.uyeq[:,:] + ini.c_eq[1]*ini.uysq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,3] = fct2*(1. - ini.c_eq[0]*ini.uxeq[:,:] + ini.c_eq[1]*ini.uxsq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,4] = fct2*(1. - ini.c_eq[0]*ini.uyeq[:,:] + ini.c_eq[1]*ini.uysq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,5] = fct3*(1. + ini.c_eq[0]*ini.uxuy5[:,:] + ini.c_eq[1]*ini.uxuy5[:,:]*ini.uxuy5[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,6] = fct3*(1. + ini.c_eq[0]*ini.uxuy6[:,:] + ini.c_eq[1]*ini.uxuy6[:,:]*ini.uxuy6[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,7] = fct3*(1. + ini.c_eq[0]*ini.uxuy7[:,:] + ini.c_eq[1]*ini.uxuy7[:,:]*ini.uxuy7[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,8] = fct3*(1. + ini.c_eq[0]*ini.uxuy8[:,:] + ini.c_eq[1]*ini.uxuy8[:,:]*ini.uxuy8[:,:]  - ini.c_eq[2]*ini.usq[:,:])

    #Collision step

    ini.tau[:,:,:] = maximum(ini.tau0, (1 - (ini.feq[:,:,:]/ini.ftemp[:,:,:])))
    for a in range(9):
       ini.f[1:ini.sizeX_+1,:,a] = np.where(ini.m[1:ini.sizeX_+1,:] == 0, ini.ftemp[1:ini.sizeX_+1,:,a] - (ini.ftemp[1:ini.sizeX_+1,:,a] - ini.feq[1:ini.sizeX_+1,:,a]) / ini.tau[1:ini.sizeX_+1,:,a], ini.f[1:ini.sizeX_+1,:,a])

    #Plotting 
    print "Time = ", t
    print "Mass = ", sum(ini.rho)
    print "Velocity x dir = ", sum(ini.ux)
    print "Velocity y dir = ", sum(ini.uy)
    print ini.f[2:ini.sizeX_,2:ini.sizeY_,:].min()
    #To save the density and velocity distribution
    # np.save(os.path.join(ini.name, "rho_" + ini.name + "_" + str(t).zfill(4)), ini.rho)
    # np.save(os.path.join(ini.name, "ux_" + ini.name + "_" + str(t).zfill(4)), ini.ux)     

    #Plotting the heatmaps
    if mod(t,100) == 0:
        varm = ini.u.transpose()        #Change the variable to the one that will be plotted: rho, ux, or uy
        plt.figure(1)
        ax = sns.heatmap(varm, annot=False, vmin=0., vmax=0.2, cmap='RdYlBu_r')           #For ux
        #figure(num=1, figsize=(15,6), dpi=80, facecolor='w', edgecolor='k')
        fig = plt.gcf()
        fig.set_size_inches(20.,5., forward=True)
        #ax = sns.heatmap(varm, annot=False, vmin=5.*ini.f_init, vmax=18*ini.f_init, cmap='RdYlBu_r')           #For rho
        ax.invert_yaxis()
        #plt.pause(0.001)
        plt.savefig("vel."+ini.name2+"_"+str(t).zfill(6)+".png")
        plt.clf()
    
    #Plot drag force vs time
    # if mod(t,5) == 0:
    #     plt.figure(3)
    #     plt.plot(range(ini.T), ini.Ft, label='Outlet')
    #     plt.legend(loc='upper right')
    #     plt.show
    #     plt.pause(0.001)
    #     plt.clf()

#To save the drag force
# np.save(os.path.join(ini.name, "dragForce_" + ini.name), ini.Ft)

####################################################### OUTPUT ###########################################################################