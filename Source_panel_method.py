# These libraries are required to run the code
# matplotlib
import matplotlib.pyplot as plt
# numpy
import numpy as np
# math
import math
# scipy
from scipy import integrate

# This is the main function.
def main():
	nums = [4,8,16,32]
	for num in nums:
		ax1, ax2 = plot_setup()
		panelMethod(num, ax1, ax2)
		ax1.set_title('Panel Method Shape, Number of Panels: '+ str(num), fontsize='24')
		plt.show()


# Defining a Panel object, as there are many and this object needs to be reused multiple times.
# Thus it is good to create a class (or object) out of it
class Panel:
	# constructor method to initialize different parameters for each Panel
	def __init__(self, x0, y0, x1, y1 ):
		self.x0 = x0
		self.y0 = y0
		self.x1 = x1
		self.y1 = y1
		self.len = math.sqrt((x1-x0)**2 + (y1-y0)**2)
		self.x_cent = (x0 + x1) / 2
		self.y_cent = (y0 + y1) / 2
		#self.beta = -(math.atan((y1 - y0)/(x1 - x0)) + math.pi)
		if x1-x0 <= 0.:
			self.beta = math.acos((y1-y0)/self.len)
		elif x1-x0 > 0.:
			self.beta = math.pi + math.acos(-(y1-y0)/self.len)
		self.U_theta = 0
		self.cp = 0
		self.sig = 0

# creates the panel plots
def panelMethod(num_panels, ax1, ax2):
	
	R = 1.0
	alop, x_cor, y_cor = createPanels(num_panels, R)
	theta = np.linspace(0, 2*math.pi, 100)
	x_cyl, y_cyl = R*np.cos(theta), R*np.sin(theta)
	plot_shape(alop, x_cor, y_cor, ax1)
	calcCp(ax2, num_panels, alop, 1.0, R)
	plotCp(ax2, num_panels, alop, 1.0, x_cyl, y_cyl, R)

# creates a symmetric polygon from a specified number of panels
def createPanels(num_panels, R):
	# define the end-points of the panels
	x_cor = R*np.cos(np.linspace(0, 2*math.pi, num_panels+1))
	y_cor = R*np.sin(np.linspace(0, 2*math.pi, num_panels+1))

	# define the panels
	alop = np.empty(num_panels, dtype=object)
	for i in range(num_panels):
		alop[i] = Panel(x_cor[i], y_cor[i], x_cor[i+1], y_cor[i+1])

	return alop, x_cor, y_cor

# compute the integral of the coefficient of pressure components, as in the book (Eq. 3.159)
def integralCp(CpX, CpY):
	def integr(x):
		# Equation 3.159 in the book
		eq1 = (CpX.y_cent - (math.cos(CpY.beta) * x + CpY.y1)) * math.cos(CpX.beta)
		eq2 = (-(CpX.x_cent - (-math.sin(CpY.beta) * x + CpY.x1)) * math.sin(CpX.beta))
		eq3 = CpX.x_cent - (CpY.x1 - math.sin(CpY.beta) * x)
		eq4 = CpX.y_cent - (CpY.y1 + math.cos(CpY.beta) * x)
		return (eq1 + eq2) / (eq3**2 + eq4**2)
	return integrate.quad(integr, 0, CpY.len)[0]

# Here, all the calculations to set the panel parameters using lin. alg.
def calcCp(ax, num_panels, alop, U, R):
	# zero matrix
	matrixToSolve = np.zeros((num_panels, num_panels))

	# adding computed (integrated) values to all elements in the above matrix, except for the diagonal
	for i in range(0,len(alop)):
		for j in range(0,len(alop)):
			if i != j:
				CpX = alop[i]
				CpY = alop[j]
				matrixToSolve[i,j] = 0.5 / math.pi * integralCp(CpX, CpY)

	b = - U * np.sin([pan.beta for pan in alop])
	sig = np.linalg.solve(matrixToSolve, b)

	for i in range(0,len(alop)):
		alop[i].sig = sig[i]

	U_theta = np.dot(matrixToSolve, sig) + b

	for i in range(0,len(alop)):
		alop[i].U_theta = U_theta[i]
		alop[i].cp = 1.0 - (alop[i].U_theta/U)**2

# plot setup
def plot_setup():
	fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(4,4), facecolor='white')
	ax2.set_title('Analytical vs. Panel Method', fontsize=24)
	ax1.set_xlabel('x', fontsize='16'), ax1.set_ylabel('y', fontsize='16')
	ax1.grid()
	ax2.grid()
	return ax1, ax2

# plot the shape (of panels)
def plot_shape(alop, x_cor, y_cor, ax):
	ax.plot(x_cor, y_cor, linestyle='-')
	ax.scatter([p.x0 for p in alop], [p.y0 for p in alop])
	ax.scatter([p.x_cent for p in alop], [p.y_cent for p in alop], color='k')

# create the Cp vs. X plot
def plotCp(ax, num_panels, alop, u_inf, x_cylinder, y_cylinder, R):
	#plotting
	# calculate the analytical surface pressure coefficient
	cp_analytical = 1.0 - 4*(y_cylinder/R)**2

	# plot the surface pressure coefficient
	ax.set_xlabel('x', fontsize='16')
	ax.set_ylabel('$C_p$', fontsize='16')
	ax.plot(x_cylinder, cp_analytical, color='b')
	ax.scatter([p.x_cent for p in alop], [p.cp for p in alop])
	ax.set_xlim(-1.1, 1.1)
	ax.set_ylim(-4, 2);

# run the main function
main()
