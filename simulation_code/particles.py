from numpy import *
from pylab import rand
import pandas
import numpy as np

import matplotlib.pyplot as plt
reload(plt)
reload(pandas)

sturges = lambda n: int(log2(n) + 1)
square_root = lambda n: int(sqrt(n))
from scipy.stats import kurtosis
doanes = lambda data: int(1 + log(len(data)) + log(1 + kurtosis(data) * (len(data) / 6.) ** 0.5))

# PARTICLES MODULE:
# This set of classes was developed to model interacting spherical particles in 3
# dimensions. There are presently two types of interaction forces possible. Each
# is a class
#
# GranularMaterialsForce: This is a force that represents elastic forces between
# particles that are overlapping, and a damping, or dissipative force. It does
# not include the often used rotational degrees of freedom. Also note that the
# constraint forces are taken care of with an extra function call.
#
# LennardJonesForce: This force describes the weak attraction at large distances
# and strong repulsion experienced at close distances by mono-atomic gases.
#
# The equation of motion are integrated by the Verlet method, which is presently
# the only integration scheme supported.
#
# All particle data (position, velocity, acceleration, radii, and distances) are
# stored in the Particles class.

# 9/23/08 JVJ and Tim Bocek

class GranularMaterialForce:

    def __init__(self, k = 1.5, gamma = .3, g = .55):
        # parameters in force model
        self.__k = k          # Elastic 'bounce'
        self.__gamma = gamma  # Energy dissipation/loss
        self.__g = g          # Gravity

    def __call__(self, p):
        # Find position differences
        d, dx, dy, dz = p.distanceMatrix(p.x,p.y,p.z)

        # Compute overlap
        dr = d - triu(p.sumOfRadii) - tril(p.sumOfRadii)

        # No forces arising in no overlap cases
        dr[dr>0]=0

        # Compute elastic particle/particle forces
        magnitude = -self.__k * dr

        # Velocity differences
        dv, dvx, dvy, dvz = p.distanceMatrix(p.vx,p.vy,p.vz)

        # Damping terms
        vijDotrij = dvx * dx + dvy * dy + dvz * dz
        vijDotrij[dr==0] = 0.

        # Damping is subtracted from force
        magnitude = magnitude - self.__gamma *  vijDotrij / d

        magnitude[d==-1] = 0

        cx, cy, cz = self.floorConstraint(p)

        # Project onto components, sum all forces on each particle
        p.ax =  sum(magnitude * (-dx/d) * p.ratioOfRadii,axis=1) + cx
        p.ay =  sum(magnitude * (-dy/d)* p.ratioOfRadii,axis=1) - self.__g + cy
        #p.netSum = (dy[dy<1])
        #print("netSum")
        #print(p.netSum)
        p.az =  sum(magnitude * (-dz/d)* p.ratioOfRadii,axis=1) + cz

    def floorConstraint(self,p):
        " "" This is a highly specific function for a floor that responds (elasticity and damping) the same way a particle does. Presently, if constraints are to change, one would have to rewrite the function."" "
        effectiveRadius = 3. # This is how 'hard' the floor is
        floorDistance = p.y + p.L/2 - p.r
        floorDistance[floorDistance > 0] = 0
        lowerWallForce = -self.__k * floorDistance
        lowerWallDamping = -self.__gamma * p.vy * floorDistance
        lowerWallForce = lowerWallForce - lowerWallDamping
        cx = 0
        cy = lowerWallForce * effectiveRadius / p.r
        cz = 0
        return cx, cy, cz

class VerletIntegrator:

    def __init__(self,dt=0.01):
        # Time step
        self.__dt = dt

    def __call__(self,force,p):


        # Position update
        p.x = p.x + p.vx * self.__dt + .5 * p.ax * self.__dt ** 2
        p.y = p.y + p.vy * self.__dt + .5 * p.ay * self.__dt ** 2
        p.z = p.z + p.vz * self.__dt + .5 * p.az * self.__dt ** 2

        # Update periodic BC
        p.pbcUpdate()

        # Store accelerations for averaging that is done
        ax=p.ax
        ay=p.ay
        az=p.az

        force(p) # Force update with new positions

        # Velocity updates
        p.vx = p.vx + 0.5 * (ax + p.ax) * self.__dt
        p.vy = p.vy + 0.5 * (ay + p.ay) * self.__dt
        p.vz = p.vz + 0.5 * (az + p.az) * self.__dt


class Particles:

    def __init__(self,L,force,periodicX = 1,periodicY=1,periodicZ=1):
        # Container size
        self.L = L
        # Total Number of particles
        self.N = 0
        # type
        self.type = 'float32'
        # Positions
        self.x = array([],dtype=self.type)
        self.y = array([],dtype=self.type)
        self.z = array([],dtype=self.type)
        # Velocities
        self.vx = array([],dtype=self.type)
        self.vy = array([],dtype=self.type)
        self.vz = array([],dtype=self.type)
        # Forces
        self.ax = array([],dtype=self.type)
        self.ay = array([],dtype=self.type)
        self.az = array([],dtype=self.type)
        # Radii
        self.r = array([],dtype=self.type)
        # Periodic on?
        self.periodicX = periodicX
        self.periodicY = periodicY
        self.periodicZ = periodicZ
        # Force function
        self.f = force

    def log(self):
        print ("Radiuses:")
        print(self.r)

        print("Acceleration in y")
        print(self.ay)

        print("Y position")
        print(self.y)

        print("Velocity in y")
        print(self.vy)

        d = {'z': np.around(self.z, decimals=5), 'x': np.around(self.x, decimals=4), 'radius': np.around(self.r, decimals=4), 'ay': np.around(self.ay, decimals=4), 'y': np.around(self.y, decimals=4), 'vy': np.around(self.vy, decimals=4)}
        nutz = pandas.DataFrame(d)
        print(nutz)
        nutz.to_csv("nutz.csv", index=False, cols=['radius', 'ay', 'z', 'y', 'x', 'vy'])
        #nutz.y.hist()

        #bp = nutz.boxplot(column="y", by="radius", grid=False)
        #plt.ylabel("[y position]")
        #for i in [1,3]:
        #    y = nutz.y[nutz.radius==i].dropna()
        #    # Add some random "jitter" to the x-axis
        #    e = i
        #    if (i==3):
        #        e=2
        #    x = np.random.normal(e, 0.04, size=len(y))
        #    print str(i) + ":" + str(x)
        #    print y
        #    plt.plot(x, y, 'r.', alpha=0.2)
        #    #plt.show()
        #plt.savefig("boxplot.png", dpi=300)

        #bp2 = nutz.boxplot(whis=5, column="ay", by="radius", grid=False)
        #plt.ylabel("[Sum of forces in y direction]")
        #for i in [1,3]:
        #    y = nutz.ay[nutz.radius==i].dropna()
        #    e = i
        #    if (i==3):
        #        e=2
        #    x = np.random.normal(e, 0.04, size=len(y))
        #    plt.plot(x, y, 'r.', alpha=0.2)
        #plt.savefig("boxplot2.png", dpi=300)

        #titanic.fare.hist(bins=doanes(titanic.fare.dropna()), normed=True, color='lightseagreen')
        plt.figure()
        plt.xlabel('Sum of forces in y direction')
        hist1 = nutz.ay.hist(bins=doanes(nutz.ay.dropna()), normed=True, grid=False, color='lightseagreen')
        #hist1 = nutz.ay.hist(bins=30, normed=True, grid=False, color='lightseagreen')
        #titanic.fare.dropna().plot(kind='kde', xlim=(0,600), style='r--')
        nutz.ay.dropna().plot(kind='kde', style='r--')
        plt.savefig("hist_doanes.png", dpi=300)

        plt.figure()
        plt.xlabel("Sum of forces in y direction")
        hist2 = nutz.ay.hist(bins=30, normed=True, grid=False, color='lightseagreen')
        nutz.ay.dropna().plot(kind='kde', style='r--')
        plt.savefig('hist_30.png', dpi=300)

    def addParticle(self,x,y,z,vx,vy,vz,r):
        self.x = hstack((self.x,x))
        self.y = hstack((self.y,y))
        self.z = hstack((self.z,z))
        self.vx = hstack((self.vx,vx))
        self.vy = hstack((self.vy,vy))
        self.vz = hstack((self.vz,vz))
        self.r = hstack((self.r,r))
        self.N = self.N+1
        temp = tile(self.r,(self.N,1))
        self.sumOfRadii = temp + temp.T
        self.ratioOfRadii = temp / temp.T
        self.f(self)

    def pbcUpdate(self):
        " ""Moves paricles across periodic boundary"" "
        if self.periodicX:
            self.x[self.x > self.L/2]  = self.x[self.x > self.L/2]  - self.L
            self.x[self.x < -self.L/2] = self.x[self.x < -self.L/2] + self.L
        if self.periodicY:
            self.y[self.y > self.L/2]  = self.y[self.y > self.L/2]  - self.L
            self.y[self.y < -self.L/2] = self.y[self.y < -self.L/2] + self.L
        if self.periodicZ:
            self.z[self.z > self.L/2]  = self.z[self.z > self.L/2]  - self.L
            self.z[self.z < -self.L/2] = self.z[self.z < -self.L/2] + self.L

    def distanceMatrix(self,x,y,z):
        " ""Computes distances between all particles and places the result in a matrix such that the ij th matrix entry corresponds to the distance between particle i and j"" "
        xtemp = tile(x,(self.N,1))
        dx = xtemp - xtemp.T
        ytemp = tile(y,(self.N,1))
        dy = ytemp - ytemp.T
        ztemp = tile(z,(self.N,1))
        dz = ztemp - ztemp.T

        # Particles 'feel' each other across the periodic boundaries
        if self.periodicX:
            dx[dx>self.L/2]=dx[dx > self.L/2]-self.L
            dx[dx<-self.L/2]=dx[dx < -self.L/2]+self.L
        if self.periodicY:
            dy[dy>self.L/2]=dy[dy>self.L/2]-self.L
            dy[dy<-self.L/2]=dy[dy<-self.L/2]+self.L
        if self.periodicZ:
            dz[dz>self.L/2]=dz[dz>self.L/2]-self.L
            dz[dz<-self.L/2]=dz[dz<-self.L/2]+self.L

        # Total Distances
        d = sqrt(dx**2+dy**2+dz**2)

        # Mark zero entries with negative 1 to avoid divergences
        d[d==0] = -1


        #print "Dist matrix dy"
        #print dy
        return d, dx, dy, dz
