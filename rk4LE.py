#!/usr/bin/python
# -*- coding: latin-1 -*-

#This script solves the Lane-Emden equation numerically. The equation is a toy model of a star.
#The solution is returned in terms of ascii output with headers. First column is radius in meters.
#Second column is density in kg/m^3. Third column is pressure in Pascal. Internal (dimensionless)
#quantities are used in the computation. Fourth column is dimensionless radius, usually known as xi.
#Fifth colum is the dimensionless quantity theta. Sixt column its derivative wrt dimensionless radius.
#Run using ./rk4_laneemden.py > lane. The file "lane" can be read using read.table() in R, for example.
#the script will run in a linux or mac-os system with python installed in /usr/bin/python.
#If you cannot run it check that it is executable. chmod u+x rk4_laneemden.py will make it executable.
#Read more about the equation at http://en.wikipedia.org/wiki/Lane–Emden_equation
#This script uses the fourth-order Runge-Kutta algorithm to solve the equation numerically.
#I know that there are better algorithms around, this was written to show RK4 to computational astronomy
#students. It may not be optimal but it was written together with the students in an application
#of the socratic method of teaching (maieutics) and was only minimally edited afterwards.

#Thinking points:
#Solving an ordinary differential equation is possible if we provide initial conditions.
#This is true both for analytic solutions and for numeric solutions. This is a second order equation, so
#we need to provide both an initial value for the function and an initial value for its derivative.
#The equation does not "know" anything about the physics of the problem. It may happen that the dimensionless
#variable theta becomes negative for some radius, which would mean that the density is negative. This
#is physically unacceptable but the program does not know it. We have to realize it for ourselves.
#If we want our "star" to be similar to the sun we can set the proportionality constant K in the polytropic
#equation of state and the central density so that both density and pressure are the same as those in the
#center of the sun. But how do we know these, except from Wolfram alpha or a book? Shouldn't we use the
#surface density and pressure instead? The most straightforward way to provide initial conditions for an
#equation may not be the best from the operational point of view!
#What about the equation of state of ideal gases? Why can we use this polytropic equation instead?
#Is the Sun a perfect gas? (spoiler: no, it has a different equation of state because of radiation)
#Does this equation describe the correct physics of energy generation in the sun? (no...)
#If not, why does it (kind-of) work?

import math

n = 5.0 #Lane-Emden exponent. Analytical solutions exist for n = 0, 1, 5.
R = 60.0 #Dimensionless radius where to stop the computation. Increase to get an error.
gammapoly = (1.0 + (1.0/n))

#Dimensional constants are not used in the computation. We use them just to scale the output.
#Why is this a good idea? Why is computing with big numbers a bad idea? Does this machine precision stuff
#apply to Python as well? Figure this stuff out as an exercise...
K = 5000000000.0 #pressure P = K rho^(1 + 1/n). 
rho_c = 100000.0 #central density kg/m^3.
G = 6.67*math.pow(10.0, -11.0) #gravitational constant
alpha = math.pow(K*(n+1)*math.pow(rho_c, gammapoly - 2.0)/(4*math.pi*G), 0.5) #derived scale length

#This function returns the second derivative of theta as a function of xi.
#v = dtheta/dxi, x = theta, t = xi
#this misleading and useless notation was left here to show the parallel with mechanics:
#we have something similar to evolution in a potential plus a friction force that decreases over time
def a(v, x, t):
	return - math.pow(x, n) - 2.0*v/t

#initial conditions
x = 1.0 #theta is 1 in the center. Changing the value changes the scale for density.
rho = rho_c #dimensional density in the center is equal to the dimensional density in the center
P = K*math.pow(rho, gammapoly) #dimensional pressure in the center obeys the polytropic equation of state
v = 0.0 #dtheta/dxi is 0 in the center. We do not have a cusp in density in the center.
h = 0.001 #stepsize in dimensionless radius. The smaller the more accurate and the slower.
t = 0.0001 #initial dimensionless radius. Why not 0.0?
r = alpha*t #dimensional radius

print "r rho P xi theta dthetadxi" #just an header for making it easier (or more difficult?) to read the file
while t <= R: #the RK loop, finally
	print r, rho, P, t, x, v #we print the dimensional quantities first (radius, density, pressure) then the others
	#the RUNGE-KUTTA thing!
 	k1_x = v #estimate of the velocity
	k1_v = a(v, x, t) #let's use this to estimate acceleration
	k2_x = v + 0.5*h*k1_v #so we get a different estimate of the velocity
	k2_v = a(v + 0.5*h*k1_v, x + 0.5*h*k1_x, t + 0.5*h) #another estimate of acceleration
	k3_x = v + 0.5*h*k2_v #another estimate of velocity
	k3_v = a(v + 0.5*h*k2_v, x + 0.5*h*k2_x, t + 0.5*h) #you get the idea
	k4_x = v + h*k3_v #...
	k4_v = a(v + h*k3_v, x + h*k3_x, t + h) #...
	x = x + h*(k1_x + 2.0*k2_x + 2.0*k3_x + k4_x)/6.0 #putting them all together
	#RK magic is over
	rho = rho_c*math.pow(x, n) #dimensional density
	r = alpha*t #dimensional radius
	P = K*math.pow(rho, gammapoly) #dimensional pressure
	#update velocity and time
	v = v + h*(k1_v + 2.0*k2_v + 2.0*k3_v + k4_v)/6.0
	t = t + h

#We are done!

