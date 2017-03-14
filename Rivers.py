import random as random
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline

# Insert the initial two endpoints, at (0,0) and (1,0).
x1 = 0
y1 = 0
x2 = 1
y2 = 0p1 = [x1,y1]
p2 = [x2,y2]
"
# Create the list of points, and add the first two points to it.
P = []
P.append(p1)
P.append(p2)

temp = [a for a in P]

"N = 1 # The number of initialy subintervals. This code currently doesn't
      # have functionality for more than one built in.

def distance(point1,point2):
    """
    Determine the distance between two points.
    
    This function is used later on when parameterizing the curve.
    
    Parameters:
    -----------
    point1 : list
        first point
    point2 : list
        second point
    
    Returns:
    --------
    float
        Distance between the two points
    """
    x1 = point1[0]
    y1 = point1[1]
    x2 = point2[1]
    y2 = point2[1]
    xdiff = x1 - x2
    ydiff = y1 - y2
    dist = np.sqrt(xdiff**2 + ydiff**2)
    return dist

for j in range(0,4):
    def add_point(i):
        """
        Create a new point and add it to a temporary list of points.

        Parameters:
        -----------
        i : int
            The number of the iteration minus one
        Returns:
        --------
        None
        """
        p1 = P[i]
        p2 = P[i+1]
        x1 = p1[0]
        y1 = p1[1]
        x2 = p2[0]
        y2 = p2[1]

        xdiff = x2 - x1
        ydiff = y2 - y1
        m = ydiff/xdiff # Slope of line connecting the two points
        b = y1 - m*x1
        xmid = x1 + xdiff/2
        ymid = y1 + ydiff/2
        midpoint = [xmid,ymid] # Midpoint of line segment
        if m != 0:
            mnew = -m**(-1) # Slope of line perpendicular to the original line segment
            theta = np.arctan(mnew)
        else:
            theta = np.pi/2
        r = random.uniform(-.1/(j+1),.1/(j+1))
        xnew = xmid + r*np.cos(theta)
        ynew = ymid + r*np.sin(theta)
        pnew = [xnew,ynew] # New point
        temp.insert(i*2+1,pnew)

    for i in range(0,len(P) - 1):
        add_point(i)

    P = [a for a in temp]

X = []
Y = []
for i in range(0,len(P)):
    point = P[i]
    x = point[0]
    y = point[1]
    X.append(x)
    Y.append(y)

T = [0]

for i in range(0,len(P)-1):
    dist = distance(P[i],P[i+1])
    T.append(dist) # Parameterizes the curve

Tnew = np.linspace(T[0],T[len(T)-1],100000)
Xsmooth = spline(T,X,Tnew)
Ysmooth = spline(T,Y,Tnew)

plt.plot(Xsmooth,Ysmooth)
plt.xlim(0,1)
#plt.plot(X,Y) # Optional
plt.show()
