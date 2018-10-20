import matplotlib.pyplot as plt
import numpy as np 
import os

#Plotting dragforce in the cylinder

q = np.arange(4)
print q[:]

T = 200
t = range(1000)
x = np.load(os.path.join("poiseuille33_ux001", "ux_poiseuille33_ux001_0190.npy"))
y = np.load(os.path.join("poiseuille33_ux003", "ux_poiseuille33_ux003_0190.npy"))
a = np.load(os.path.join("poiseuille33_ux005", "ux_poiseuille33_ux005_0190.npy"))
b = np.load(os.path.join("poiseuille33_ux01", "ux_poiseuille33_ux01_0190.npy"))
c = np.load(os.path.join("poiseuille33_ux02", "ux_poiseuille33_ux02_0190.npy"))
d = np.load(os.path.join("poiseuille33_ux03", "ux_poiseuille33_ux03_0190.npy"))
e = np.load(os.path.join("poiseuille33_ux04", "ux_poiseuille33_ux04_0190.npy"))
#d = np.load(os.path.join("r45_u3", "dragForce_r45_u3.npy"))
sizeX_ = len(a) - 2         #length in x-direction
sizeY_ = len(a[0]) - 2        #length in y-direction

plt.figure(1)
plt.plot(x[(sizeX_+2)//2,:], label='u1 = 0.01')
plt.plot(y[(sizeX_+2)//2,:], label='u2 = 0.03')
plt.plot(a[(sizeX_+2)//2,:], label='u3 = 0.05')
plt.plot(b[(sizeX_+2)//2,:], label='u4 = 0.1')
plt.plot(c[(sizeX_+2)//2,:], label='u5 = 0.2')
plt.plot(d[(sizeX_+2)//2,:], label='u6 = 0.3')
plt.plot(e[(sizeX_+2)//2,:], label='u7 = 0.4')
#plt.plot(t, d, label='r4')
plt.legend(loc='upper right')
plt.savefig('u_vs_y.eps', format='eps', dpi=1000000)
plt.show()
