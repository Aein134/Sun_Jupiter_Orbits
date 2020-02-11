import numpy as np
import matplotlib.pyplot as plt
import sys

###Call with python ./HW2.py x_arg y_arg u_arg v_arg case#


points = 100000
h = 0.01

x0 = sys.argv[1]
y0 = sys.argv[2]
u0 = sys.argv[3]
v0 = sys.argv[4]
case = sys.argv[5]


def xdot(x, y, u, v):
 return u

def ydot(x, y, u, v):
 return v

def udot(x, y, u, v):
 Msun = 1.989*(10**30)
 Mjupiter = 1.898*(10**27)
 MuSun = Msun/(Msun+Mjupiter)
 MuJupiter = Mjupiter/(Msun+Mjupiter)
 rSun = ((x - MuJupiter)**2 + y**2)**(1.0/2.0)
 rJupiter = ((x + MuSun)**2 + y**2)**(1.0/2.0)
 return -MuSun * (x - MuJupiter)/(rSun**3) - MuJupiter * (x + MuSun)/(rJupiter**3) + x + 2*v

def vdot(x, y, u, v):
 Msun = 1.989*(10**30)
 Mjupiter = 1.898*(10**27)
 MuSun = Msun/(Msun+Mjupiter)
 MuJupiter = Mjupiter/(Msun+Mjupiter)
 rSun = ((x - MuJupiter)**2 + y**2)**(1.0/2.0)
 rJupiter = ((x + MuSun)**2 + y**2)**(1.0/2.0)
 return -MuSun * (y)/(rSun**3) - MuJupiter * (y)/(rJupiter**3) + y - 2*u

def F(x, y, u, v, n):
 if n == 1:
  return xdot(x, y, u, v)
 elif n == 2:
  return ydot(x, y, u, v)
 elif n == 3:
  return udot(x, y, u, v)
 elif n == 4:
  return vdot(x, y, u, v)

K = np.zeros([4,4])

x = np.empty([points])
y = np.empty([points])
u = np.empty([points])
v = np.empty([points])
x[0] = x0
y[0] = y0
u[0] = u0
v[0] = v0

for i in range(points-1):
#################################### Calculate the K values #####################################
 for j in range(4):
  for k in range(4):
   if j+1 == 1:
    K[j,k] = h*F(x[i], y[i], u[i], v[i], k+1)
   elif j+1 == 2:
    K[j,k] = h*F(x[i] + K[0,0]/2.0, y[i] + K[0,1]/2.0, u[i] + K[0,2]/2.0, v[i] + K[0,3]/2.0, k+1)
   elif j+1 == 3:
    K[j,k] = h*F(x[i] + K[1,0]/2.0, y[i] + K[1,1]/2.0, u[i] + K[1,2]/2.0, v[i] + K[1,3]/2.0, k+1)
   elif j+1 == 4:
    K[j,k] = h*F(x[i] + K[2,0], y[i] + K[2,1], u[i] + K[2,2], v[i] + K[2,3], k+1)
#################################################################################################

 x[i+1] = x[i] + (1.0/6.0)*(K[0,0] + 2*K[1,0] + 2*K[2,0] + K[3,0])
 y[i+1] = y[i] + (1.0/6.0)*(K[0,1] + 2*K[1,1] + 2*K[2,1] + K[3,1])
 u[i+1] = u[i] + (1.0/6.0)*(K[0,2] + 2*K[1,2] + 2*K[2,2] + K[3,2])
 v[i+1] = v[i] + (1.0/6.0)*(K[0,3] + 2*K[1,3] + 2*K[2,3] + K[3,3])

 if i%10000 == 0:
  print(i)

plt.plot(x, y)
plt.savefig("jupiter_sun_orbit_case"+case+".png")
plt.show()
