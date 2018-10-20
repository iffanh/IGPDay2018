####################################################### PREPARATION ###########################################################################
from numpy import *
import scipy.io as sp 


# #importing grid
m = genfromtxt('../porestructure/channel_33.dat', delimiter="\t")
m = m.transpose()

#for a new file, create folder named after this variable name
name = "poiseuille33_ux003"

#What we have now:
# 1. channel_5.dat
# 2. channel_9.dat
# 3. channel_33.dat
# 4. cylinder_5.dat
# 5. cylinder_11.dat
# 6. cylinder_15.dat

print m

sizeX_ = len(m) - 2         #length in x-direction
sizeY_ = len(m[0]) - 2        #length in y-direction

#The number of iteration
T = 200            #Total time used in the simulation
dt = 1             #time interval

#Declaring variables

rho = zeros((sizeX_+2,sizeY_+2))                         #Density of the lattice point, 
ux = zeros((sizeX_+2,sizeY_+2))                          #Macroscopic velocity of the lattice point 
uy = zeros((sizeX_+2,sizeY_+2))
u = zeros((sizeX_+2,sizeY_+2))
uxeq = zeros((sizeX_+2,sizeY_+2))                        #Macroscopic velocity of the lattice point 
uyeq = zeros((sizeX_+2,sizeY_+2))
f = zeros((sizeX_+2,sizeY_+2,9))      #Density distribution of the a point f[x position][y position][index]
ftemp = zeros((sizeX_+2,sizeY_+2,9))
feq = zeros((sizeX_+2,sizeY_+2,9))
tau = zeros((sizeX_+2,sizeY_+2,9))

#Constants used
Re_x = 3.
ux0 = 0.03                       #for constant velocity boundaries
rhoi = 2.0
rhoo = 0.1                      #for constant pressure boundaries
uy0 = 0.2   
visc = ux0*11./Re_x             #kinematic viscosity
tau0 = 3*visc + 0.5                
e_ = array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],[0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])         
w = array([4.0/9.0, 1.0/9.0, 1.0/36.0])
c_eq = array([3., 9./2., 3./2.])
f_tol = 0.01 

####################################################### FUNCTIONS ###########################################################################
#To find out whether a m is in contact in fluid or not
def is_interior(i,j):       
    if m[i,j] == 1 and m[i-1,j] == 1 and m[i+1,j] == 1 and m[i,j-1] == 1 and m[i,j+1] == 1 and m[i-1,j-1] == 1 and m[i+1,j-1] == 1 and m[i-1,j+1] == 1 and m[i+1,j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 0.10
for j in range(1,sizeY_+ 1):
    for i in range(1, sizeX_+ 1):                 
        if m[i,j] == 0:
            for a in range(9):
                f[i,j,a] = f_init 

###Von Neumann Boundary condition
#Initializing flux boundary density distribution
for i in range(1,sizeX_+ 1):
    for j in range(1,sizeY_ + 1):
        if m[i,j] == 2:
            #West side
            for a in [0,2,3,4,6,7]:
                f[i,j,a] = f_init

        # if m[i,j] == 4:
        #     #East side
        #     for a in [0,1,2,4,5,8]:
        #         f[i,j,a] = f_init

#print f[:,:,[0,2,3,4,6,7]]
###Reynold Number

print "kinematic viscosity: ", visc
print "relaxation parameter, tau: ", tau0


