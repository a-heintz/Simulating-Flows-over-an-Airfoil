# These libraries are required to run the code
# matplotlib
import matplotlib.pyplot as plt
# numpy
import numpy as np
# sympy
from sympy.abc import x,y
import sympy
from sympy import log as ln
# math
import math
# scipy
from scipy import integrate

# This is the main function.
def main():
	# set up the plots
	(ax1, ax2, ax3, ax4, ax5, ax6) = plot_setup()
	# initialize initial variables
	U = 1
	R = 2
	xo = 1
	yo = 2
	alpha = 0.0872665 # 5 degrees in radians
	# conversion for use of polar in cartesian
	theta = sympy.atan2(y-yo,x-xo) - alpha
	r = sympy.sqrt((x-xo)**2 + (y-yo)**2)
	
	# no vortex stream functions
	psi, phi = get_streamFunct(U, R, r, theta), get_potentialFunct(U,R,r,theta)
	# differentiating psi
	psiX, psiY = psi.diff(y), -psi.diff(x)
	# differentiating phi
	phiX, phiY = phi.diff(y), -phi.diff(x)
	
	# Coefficient of pressure
	theta_nv = np.linspace(0, 2*math.pi, 100)
	Cp_nv_adv = get_Cp_adv(U, theta_nv, 0) # advanced

	# plot stream function
	plot_funct(psiX,psiY, ax1, 1)
	# plot potential function
	plot_funct(phiX, phiY, ax2, 0.0001)
	# plot Coefficient of pressure
	ax3.plot(theta_nv, Cp_nv_adv, linestyle='-') # plotting advanced
	ax3.set_xlabel('$\\theta$')
	ax3.set_ylabel('$C_p$')

	# With vortex
	# Make the angle of attack 0.
	theta += alpha
	# a list of various Gamma values
	Gamma = [4*math.pi*U*R - 10,4*math.pi*U*R,4*math.pi*U*R + 10]
	# a list of subplots
	axs = [ax4,ax5,ax6]
	# run through each Gamma value
	for i in range(0,len(Gamma)):
			# create a stream function
			psi = get_streamFunct(U, R, r, theta) + Gamma[i]/(2*math.pi) * ln(r/R)
			# differentiate the stream function to get the velocity components
			psiX = psi.diff(y)
			psiY = -psi.diff(x)
			# plot stream function
			plot_funct(psiX,psiY, axs[i], 1)
			# get the coefficient of lift (with rotation)
			Cl = get_Cl(Gamma[i], R, U, 0)
			# test case to check the theoretical answer is the same as the computed answer
			#(to a certain degree of error, as the python integrator probably has some error)
			print -0.000005 < Gamma[i]/(R*U) - Cl < 0.000005

	#display the plot
	plt.show()

	# Coefficient of lift and drag test cases for non-rotating cylinder
	print -0.000005 < get_Cl(0, R, U, alpha) < 0.000005
	print -0.000005 < get_Cd(0, R, U, alpha) < 0.000005

# integrate the top half of the cylinder in the x-dir
def integralTopX(x1, x2, Gamma, R,U, alpha):
	def integr(x):
		y = math.sqrt(R**2-x**2)
		theta = math.atan2(y,x) - alpha
		u_theta = -2*U*math.sin(theta) - Gamma/(2*math.pi*R)
		return 1.0 - (u_theta/U)**2
	return integrate.quad(integr, x1, x2)[0]

# integrate the bottom half of the cylinder in the x-dir
def integralBottomX(x1, x2, Gamma, R,U, alpha):
	def integr(x):
		y = -math.sqrt(R**2-x**2)
		theta = math.atan2(y,x) - alpha
		u_theta = -2*U*math.sin(theta) - Gamma/(2*math.pi*R)
		return 1.0 - (u_theta/U)**2
	return integrate.quad(integr, x1, x2)[0]

# integrate the top half of the cylinder in the y-dir
def integralTopY(y1, y2, Gamma, R,U, alpha):
	def integr(y):
		x = math.sqrt(R**2-y**2)
		theta = math.atan2(y,x) - alpha
		u_theta = -2*U*math.sin(theta) - Gamma/(2*math.pi*R)
		return 1.0 - (u_theta/U)**2
	return integrate.quad(integr, y1, y2)[0]

# integrate the bottom half of the cylinder in the y-dir
def integralBottomY(y1, y2, Gamma, R,U, alpha):
	def integr(y):
		x = -math.sqrt(R**2-y**2)
		theta = math.atan2(y,x) - alpha
		u_theta = -2*U*math.sin(theta) - Gamma/(2*math.pi*R)
		return 1.0 - (u_theta/U)**2
	return integrate.quad(integr, y1, y2)[0]

# compute the coefficient of lift
def get_Cl(Gamma, R, U, alpha):
	top = integralTopX(-R,R,Gamma,R,U, alpha)
	bottom = integralBottomX(-R,R,Gamma,R,U, alpha)
	return 1.0/(2.0*float(R)) * (bottom - top)

# compute the coefficient of drag
def get_Cd(Gamma, R, U, alpha):
	top = integralTopY(-R,R,Gamma,R,U, alpha)
	bottom = integralBottomY(-R,R,Gamma,R,U, alpha)
	return 1.0/(2.0*float(R)) * (top - bottom)

# compute the stream function
def get_streamFunct(U, R, r, theta):
	return U * r * sympy.sin(theta) * (1 - R**2 / r**2)

# compute the potential function
def get_potentialFunct(U, R, r, theta):
	return U * r * sympy.cos(theta) * (1 + R**2 / r**2)

# compute the advanced version of the coeff. of pressure
def get_Cp_adv(U, theta, U_vortex_option):
	Uo = -2*U*np.sin(theta) + U_vortex_option
	return 1 - (Uo / U)**2

# this function makes functions more memory efficient and computationally faster
def get_lambdaFunct(func1):
	return sympy.lambdify((x,y), func1, 'numpy')

# plot setup
def plot_setup():
	fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(ncols=3, nrows=2, figsize=(4,4), facecolor='white')
	ax1.set_title('Velocity Stream Function', fontsize=18)
	ax2.set_title('Velocity Potential Function', fontsize=18)
	ax3.set_title('Cp Function vs Theta', fontsize=18)
	ax4.set_title('Stream Fn: Gamma < 4pi U R', fontsize=18)
	ax5.set_title('Stream Fn: Gamma = 4pi U R', fontsize=18)
	ax6.set_title('Stream Fn: Gamma > 4pi U R', fontsize=18)
	return (ax1, ax2, ax3, ax4, ax5, ax6)

# plotting functions in a vector field
def plot_funct(u,v, ax, arsi):
	xo = 1
	yo = 2
	R = 2
	interval = R * 3
	lowlimX, highlimX, lowlimY, highlimY = xo - interval, xo + interval, yo - interval, yo + interval
	u, v = get_lambdaFunct(u), get_lambdaFunct(v)
	Y, X = np.mgrid[lowlimX-1:highlimX+1:100j, lowlimY-1:highlimY+1:100j]
	# overlay the cylinder
	c = plt.Circle((xo, yo), radius=R, facecolor='none')
	ax.streamplot(X,Y,u(X,Y),v(X,Y), arrowsize = arsi)
	ax.set(aspect='equal')
	ax.axis([lowlimX, highlimX, lowlimY, highlimY])
	ax.add_patch(c)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')

# run the main function
main()