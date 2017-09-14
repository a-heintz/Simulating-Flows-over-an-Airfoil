import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import integrate



def joukowski(beta, a, c, alpha):
	# define constants
	i = np.sqrt(-1+0j)
	U_inf = 1
	rho = 1
	# meshing/plotting parameters
	D = 0.7
	dx = 0.01
	fs = 18
	dtheta = 0.01
	# calculate center of circle
	zeta0 = a - c * np.exp(-i*beta)
	# form matrix of values in computational domain
	something = int(round(2.0*D/dx+1))
	nx, ny = (something, something)
	xi1 = np.linspace(-D,D,nx)
	eta1 = np.linspace(-D,D,nx)
	xi, _ = np.meshgrid(xi1, eta1)
	eta = np.transpose(xi)
	zeta = xi + i * eta

	for j in range(0,len(xi)):
		for k in range(0,len(eta)):
			if abs(zeta[j][k]-zeta0) < c:
				zeta[j][k] = float('nan')

	# create subplots
	f, ((ax1, ax2, ax3),(ax4, ax5, ax6)) = plt.subplots(2, 3)
	# compute flow solution in computational domain
	gamma = 4*math.pi*U_inf*c*math.sin(alpha+beta)
	print 'Circulation:', gamma
	
	# F is complex, stores Psi and Phi, F = Phi + i*Psi
	F = (U_inf*np.exp(-i*alpha)*((zeta-zeta0) + (c**2*np.exp(2*i*alpha))/(zeta-zeta0)) + (i*gamma)/(2*np.pi)*np.log((zeta-zeta0)*np.exp(-i*alpha)/c))
	
	# plot the cylinder surface in computational domain
	theta = np.linspace(0,2*math.pi,100)
	circle = zeta0 + c*np.exp(i*theta)

	ax1.plot(np.real(circle), np.imag(circle))
	ax1.set_xlim(-D,D), ax1.set_ylim(-D,D)
	ax1.set_xlabel('$\\xi$')
	ax1.set_ylabel('$\\eta$')
	ax1.set_title('Comp. Domain')

	# plot streamlines in computational domain
	# (lines of constant Psi)
	levels = [-0.4, -0.2, 0.0, 0.2, 0.4, 0.6]
	cp2 = ax2.contour(np.real(zeta),np.imag(zeta),np.imag(F), levels, colors='k')
	
	ax2.contourf(np.real(zeta),np.imag(zeta),np.imag(F))
	ax2.clabel(cp2, inline=True, fontsize=10, colors='k')
	ax2.plot(np.real(circle), np.imag(circle))
	ax2.set_xlim(-D,D), ax2.set_ylim(-D,D)
	ax2.set_xlabel('$\\xi$', fontsize=16)
	ax2.set_ylabel('$\\eta$', fontsize=16)
	ax2.set_title('Comp. Domain')
	
	# Joukowski: perform the calculation here that
	# finds the mapped points in the z-plane
	zT = zeta + a**2 / zeta

	# plot the airfoil surface in physical domain
	airfoil = circle+a**2/circle
	ax3.plot(np.real(airfoil), np.imag(airfoil))
	ax3.set_xlim(-D,D), ax3.set_ylim(-D,D)
	ax3.set_xlabel('x', fontsize=16)
	ax3.set_ylabel('y', fontsize=16)
	ax3.set_title('Physical Domain')

	# plot streamlines in physical domain
	# (lines of constant psi)
	cp4 = ax4.contour(np.real(zT),np.imag(zT),np.imag(F), levels, colors='k')
	ax4.contourf(np.real(zT),np.imag(zT),np.imag(F))
	ax4.clabel(cp4, inline=True, fontsize=10, colors='k')
	ax4.plot(np.real(airfoil), np.imag(airfoil))
	ax4.set_xlim(-D,D), ax4.set_ylim(-D,D)
	ax4.set_xlabel('x', fontsize=16)
	ax4.set_ylabel('y', fontsize=16)
	ax4.set_title('Physical Domain')

	# find chord length
	xi_min = np.amin(np.real(circle))
	eta_min = np.imag(circle)[np.argmin(np.real(circle))]
	xi_max = np.amax(np.real(circle))
	eta_max = np.imag(circle)[np.argmax(np.real(circle))]
	x_1 = xi_min + a**2 / xi_min
	y_1 = eta_min + a**2 / eta_min
	x_2 = xi_max + a**2 / xi_max
	y_2 = eta_max + a**2 / eta_max
	chord_length = np.real(np.sqrt((x_2-x_1)**2 + (y_2-y_1)**2))
	print 'Chord Length:' , chord_length
	print 'Angle of Attack', alpha, 'rads'

	# plotting u and v
	# differentiate F
	dFdzeta = U_inf*np.exp(-i*alpha) - U_inf*a**2*np.exp(i*alpha)/(zeta - zeta0)**2 + i*gamma/(2*math.pi*(zeta - zeta0))
	# u and v
	velsU = np.real(dFdzeta)
	velsV = -np.imag(dFdzeta)
	# mag u and v
	surface_vels = np.sqrt(velsV**2+velsU**2)
	# find values on the shape
	values = []
	for j in range(0,len(surface_vels)):
		if np.isnan(surface_vels[j]).any():
			row = surface_vels[j]
			for k in range(0,len(row)-1):
				if np.isfinite(row[k]) and (np.isnan(row[k+1]) or np.isnan(row[k-1])):
					values.append((xi1[k],eta1[j],surface_vels[j][k]))
	values.sort(key=lambda x: x[0])
	xi3 = []
	eta3 = []
	vels_values = []
	for elmt in values:
		xi3.append(elmt[0])
		eta3.append(elmt[1])
		vels_values.append(elmt[2])

	xi3 = np.asarray(xi3)
	eta3 = np.asarray(eta3)
	vels_values = np.asarray(vels_values)

	# map values to physical domain & plot
	x3 = xi3 * (xi3**2+eta3**2+a**2)/(xi3**2+eta3**2)
	ax5.scatter(x3, vels_values)
	ax5.set_xlabel('x', fontsize=16)
	ax5.set_ylabel('Velocity', fontsize=16)
	ax5.set_title('Surface Velocities')

	# plot Cp as a function of x/c, one curve for top of airfoil,
	# another curve for bottom of airfoil

	Cp = 1 - (vels_values / U_inf)**2
	ax6.scatter(x3/c, Cp)
	ax6.set_xlabel('x/c', fontsize=16)
	ax6.set_ylabel('Coeff. of Pressure', fontsize=16)
	ax6.set_title('x/c vs. Cp')

	#plt.show()

	# Integrate Cp on the airfoil to find the Lift
	# & Compare with Theory

	Cl = get_Cl(gamma, a, U_inf, alpha)
	Lift = Cl * 0.5 * rho * U_inf**2
	AnalyticalLift = rho*U_inf*gamma
	print 'Lift:', Lift
	print 'Analytical Lift:', AnalyticalLift

	return AnalyticalLift/(0.5 * rho * U_inf**2)
	
# compute the coefficient of lift
def get_Cl(Gamma, R, U, alpha):
	top = integralTopX(-R,R,Gamma,R,U,alpha)
	bottom = integralBottomX(-R,R,Gamma,R,U, alpha)
	return 1.0/(2.0*float(R)) * (bottom - top)
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

joukowski(0.35, 0.1, 0.115, 0)
