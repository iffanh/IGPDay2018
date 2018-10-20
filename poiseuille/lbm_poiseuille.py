#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
from numpy import *
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import lbm_inputdata_poiseuille as ini             #import from theinput file
import os

####################################################### SIMULATION ###########################################################################

for t in range(ini.T):

    #Zou and He velocity BCs on west side
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i,j] == 2:
                ini.rho[i,j] = (ini.f[i,j,0] + ini.f[i,j,2] + ini.f[i,j,4] + 2.*(ini.f[i,j,3] + ini.f[i,j,7] + ini.f[i,j,6])) / (1 - ini.ux0*(1 + 1e-4*sin(2*j*pi/ini.sizeY_)))
                ru = ini.rho[i,j]*ini.ux0 
                ini.f[i,j,1] = ini.f[i,j,3] + (2./3.)*ru
                ini.f[i,j,5] = ini.f[i,j,7] + (1./6.)*ru - (1./2.)*(ini.f[i,j,2] - ini.f[i,j,4])
                ini.f[i,j,8] = ini.f[i,j,6] + (1./6.)*ru - (1./2.)*(ini.f[i,j,4] - ini.f[i,j,2])
    
    # ... computing density for imaging
    ini.rho[:,:] = 0.
    for a in range(9):
        ini.f[:,:,a] = np.where(ini.f[:,:,a] > ini.f_tol, ini.f[:,:,a], ini.f_tol)
        ini.rho[:,:] += np.where(ini.m[i,j] == 0, ini.f[:,:,a], 0)

    #Streaming step
    for j in range(1,ini.sizeY_+ 1):
        j_n = (j-1) 
        j_p = (j+1)
        for i in range(1,ini.sizeX_+ 1):
            i_n = (i-1) 
            i_p = (i+1)
            if (ini.m[i,j] <> 1):
                #For streaming part, if adjacent grid is a m, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction
                
                ini.ftemp[i,j,0] = ini.f[i,j,0]
                
                if (ini.m[i_p,j] <> 1):    ini.ftemp[i_p,j,1] = ini.f[i,j,1]
                else:                           ini.ftemp[i,j,3] = ini.f[i,j,1]

                if (ini.m[i,j_p] <> 1):    ini.ftemp[i][j_p][2] = ini.f[i][j][2]
                else:                           ini.ftemp[i][j][4] = ini.f[i][j][2]

                if (ini.m[i_n,j] <> 1):    ini.ftemp[i_n,j,3] = ini.f[i,j,3]
                else:                           ini.ftemp[i,j,1] = ini.f[i,j,3]

                if (ini.m[i,j_n] <> 1):    ini.ftemp[i,j_n,4] = ini.f[i,j,4]
                else:                           ini.ftemp[i,j,2] = ini.f[i,j,4]

                if (ini.m[i_p,j_p] <> 1):  ini.ftemp[i_p,j_p,5] = ini.f[i,j,5]
                else:                           ini.ftemp[i,j,7] = ini.f[i,j,5]

                if (ini.m[i_n,j_p] <> 1):  ini.ftemp[i_n,j_p,6] = ini.f[i,j,6]
                else:                           ini.ftemp[i,j,8] = ini.f[i,j,6]

                if (ini.m[i_n,j_n] <> 1):  ini.ftemp[i_n,j_n,7] = ini.f[i,j,7]
                else:                           ini.ftemp[i,j,5] = ini.f[i,j,7]

                if (ini.m[i_p,j_n] <> 1):  ini.ftemp[i_p,j_n,8] = ini.f[i,j,8]
                else:                           ini.ftemp[i,j,6] = ini.f[i,j,8]
    
    # ... and then computing macroscopic density and velocity for each lattice point, after shifting
    ini.rho[:,:] = 0.
    ini.ux[:,:] = 0.
    ini.uy[:,:] = 0.
    for a in range(9):
        ini.rho[1:ini.sizeX_+1,1:ini.sizeY_+1] += ini.ftemp[1:ini.sizeX_+1,1:ini.sizeY_+1,a]    
        ini.ux[1:ini.sizeX_+1,1:ini.sizeY_+1] += ini.e_[0,a]*ini.ftemp[1:ini.sizeX_+1,1:ini.sizeY_+1,a]
        ini.uy[1:ini.sizeX_+1,1:ini.sizeY_+1] += ini.e_[1,a]*ini.ftemp[1:ini.sizeX_+1,1:ini.sizeY_+1,a]
    ini.ux[1:ini.sizeX_+1,1:ini.sizeY_+1] = np.where(ini.rho[1:ini.sizeX_+1,1:ini.sizeY_+1] > ini.f_tol, ini.ux[1:ini.sizeX_+1,1:ini.sizeY_+1]/ini.rho[1:ini.sizeX_+1,1:ini.sizeY_+1], ini.ux[1:ini.sizeX_+1,1:ini.sizeY_+1]/(9*ini.f_tol))
    ini.uy[1:ini.sizeX_+1,1:ini.sizeY_+1] = np.where(ini.rho[1:ini.sizeX_+1,1:ini.sizeY_+1] > ini.f_tol, ini.uy[1:ini.sizeX_+1,1:ini.sizeY_+1]/ini.rho[1:ini.sizeX_+1,1:ini.sizeY_+1], ini.uy[1:ini.sizeX_+1,1:ini.sizeY_+1]/(9*ini.f_tol))
    ini.u[1:ini.sizeX_+1,1:ini.sizeY_+1] = sqrt((ini.ux[1:ini.sizeX_+1,1:ini.sizeY_+1]**2 + ini.uy[1:ini.sizeX_+1,1:ini.sizeY_+1]**2)/2)

    #Computing equilibrium distribution function
    for j in range(1,ini.sizeY_+ 1):
        for i in range(1,ini.sizeX_+ 1):
            if ini.m[i,j] == 0:
                fct1 = ini.w[0]*ini.rho[i,j]
                fct2 = ini.w[1]*ini.rho[i,j]
                fct3 = ini.w[2]*ini.rho[i,j]

                ini.uxeq[i,j] = ini.ux[i,j]                                       #uxeq will incorporate external forces, if any
                ini.uyeq[i,j] = ini.uy[i,j] 

                uxsq = ini.uxeq[i,j]*ini.uxeq[i,j]
                uysq = ini.uyeq[i,j]*ini.uyeq[i,j]

                uxuy5 = ini.uxeq[i,j] + ini.uyeq[i,j]
                uxuy6 = -ini.uxeq[i,j] + ini.uyeq[i,j]
                uxuy7 = -ini.uxeq[i,j] - ini.uyeq[i,j]
                uxuy8 = ini.uxeq[i,j] - ini.uyeq[i,j]

                usq = uxsq + uysq

                ini.feq[i,j,0] = fct1*(1.                                                   - ini.c_eq[2]*usq)
                ini.feq[i,j,1] = fct2*(1. + ini.c_eq[0]*ini.uxeq[i,j] + ini.c_eq[1]*uxsq    - ini.c_eq[2]*usq)
                ini.feq[i,j,2] = fct2*(1. + ini.c_eq[0]*ini.uyeq[i,j] + ini.c_eq[1]*uysq    - ini.c_eq[2]*usq)
                ini.feq[i,j,3] = fct2*(1. - ini.c_eq[0]*ini.uxeq[i,j] + ini.c_eq[1]*uxsq    - ini.c_eq[2]*usq)
                ini.feq[i,j,4] = fct2*(1. - ini.c_eq[0]*ini.uyeq[i,j] + ini.c_eq[1]*uysq    - ini.c_eq[2]*usq)
                ini.feq[i,j,5] = fct3*(1. + ini.c_eq[0]*uxuy5 + ini.c_eq[1]*uxuy5*uxuy5  - ini.c_eq[2]*usq)
                ini.feq[i,j,6] = fct3*(1. + ini.c_eq[0]*uxuy6 + ini.c_eq[1]*uxuy6*uxuy6  - ini.c_eq[2]*usq)
                ini.feq[i,j,7] = fct3*(1. + ini.c_eq[0]*uxuy7 + ini.c_eq[1]*uxuy7*uxuy7  - ini.c_eq[2]*usq)
                ini.feq[i,j,8] = fct3*(1. + ini.c_eq[0]*uxuy8 + ini.c_eq[1]*uxuy8*uxuy8  - ini.c_eq[2]*usq)

    #Collision step

    ini.tau[:,:,:] = maximum(ini.tau0, (1 - (ini.feq[:,:,:]/ini.f[:,:,:])))
    for a in range(9):
       ini.f[1:ini.sizeX_+1,1:ini.sizeY_+1,a] = np.where(ini.m[1:ini.sizeX_+1,1:ini.sizeY_+1] == 0, ini.ftemp[1:ini.sizeX_+1,1:ini.sizeY_+1,a] - (ini.ftemp[1:ini.sizeX_+1,1:ini.sizeY_+1,a] - ini.feq[1:ini.sizeX_+1,1:ini.sizeY_+1,a]) / ini.tau[1:ini.sizeX_+1,1:ini.sizeY_+1,a].max(), ini.f[1:ini.sizeX_+1,1:ini.sizeY_+1,a])


    #Plotting 
    print "Time = ", t
    print "Mass = ", sum(ini.rho)
    print "Velocity x dir = ", sum(ini.ux)
    print "Velocity y dir = ", sum(ini.uy)

    if mod(t,10) == 0:
        np.save(os.path.join(ini.name, "rho_" + ini.name + "_" + str(t).zfill(4)), ini.rho)
        np.save(os.path.join(ini.name, "ux_" + ini.name + "_" + str(t).zfill(4)), ini.ux)

    if mod(t,1) == 0:
        varm = ini.rho.transpose()        #Change the variable to the one that will be plotted: rho, ux, or uy
        plt.figure(1)
        #ax = sns.heatmap(varm, annot=False, vmin=0., vmax=3.0, cmap='RdYlBu_r')           #For ux
        ax = sns.heatmap(varm, annot=False, vmin=5.*ini.f_init, vmax=15*ini.f_init, cmap='RdYlBu_r')           #For rho
        ax.invert_yaxis()
        plt.pause(1.)
        plt.clf()

    # #Plot cross-section velocity
    # plt.figure(2)
    # #plt.plot(ini.ux[-2,:],label='Outlet')
    # plt.plot(ini.rho[(ini.sizeX_+2)//2 - 2,:],label='Upstream')
    # plt.plot(ini.ux[(ini.sizeX_+2)//2,:],label='Middle')
    # plt.plot(ini.rho[(ini.sizeX_+2)//2 + 2,:],label='Downstream')
    # #print ini.ux[(ini.sizeX_+2)//2,:]
    # #plt.plot(ini.ux[2,:],label='Inlet')
    # plt.legend(loc='upper right')
    # plt.show
    # plt.pause(0.001)
    # plt.clf()

    ##For calculation purposes
    #print "Upstream: ", ini.rho[(ini.sizeX_+2)//2 - 2,(ini.sizeY_+2)//2]
    #print "Downstream: ", ini.rho[(ini.sizeX_+2)//2 + 2,(ini.sizeY_+2)//2]
    #print "umax: ", ini.ux[(ini.sizeX_+2)//2 + 2,(ini.sizeY_+2)//2]

####################################################### OUTPUT ###########################################################################

        
