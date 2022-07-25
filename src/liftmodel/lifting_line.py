from scipy.optimize import fsolve
from scipy.integrate import trapezoid
from math import pi as PI
import numpy as np

# Note:
# gamma = circulation * 2 / ( wingspan * VelocityAtInfinity ) 		# divide by 4 to get circulation from matrix form

class LiftingLine:
	def __init__(self, n, chord_ratio, geom_angle, vel_ratio=(lambda x: 1) ):
		self._N = n
		self._INDEX = np.arange(n)
			# generate grid
		self._Y = np.linspace(-PI/2, PI/2, num=n)   # half revolution
		self._X = self._Y[0:-1] + 0.5*PI/(n-1)      # midpoints
			# map to sine
		self._Y = np.sin(self._Y)
		self._X = np.sin(self._X)
			# user defined functions, f(x)
		self.setChord(chord_ratio)                  # Ratio between element chord and wing span: chord(y) / wingspan
		self.setGeomAngle(geom_angle)               # Geometric angle of attack
		self.setVelRatio(vel_ratio)                 # Ratio of element velocity to velocity at infinity: Velocity(y) / Vinf
						
	# -- Private				
	
	def _integroDifferential(self, gamma, k):
		alpha = np.diff(gamma) / ( self._Y[k] - self._X )
		alpha /= (4*PI) * self._VEL_RATIO_X
		return np.sum(alpha)
		
	def _inducedAngles(self, gamma):
		func = lambda k: self._integroDifferential(gamma, k)
		angles = np.vectorize(func)(self._INDEX)
		# tip correction		
		angles[0] = self._ALPHA[0] 
		angles[-1] = self._ALPHA[-1] 
		return angles
		
	def _rootFunction(self, gamma):	
		alphai = self._inducedAngles(gamma)
		cl = 2*PI*( self._ALPHA - alphai )
		return ( self._VEL_RATIO_Y * self._CHORD )*cl - gamma
		
	# -- Public 
	
	def setChord(self, chord_ratio):
		self._CHORD = np.vectorize(chord_ratio)(self._Y)
		self._AREA = trapezoid(self._CHORD, self._Y)
		
	def setGeomAngle(self, geom_angle):
		self._ALPHA = np.vectorize(geom_angle)(self._Y)
		
	def setVelRatio(self, vel_ratio):
		vfunc = np.vectorize(vel_ratio)
		self._VEL_RATIO_Y = vfunc(self._Y)
		self._VEL_RATIO_X = vfunc(self._X)
		
	def solveCirculation(self):
		gamma = fsolve( self._rootFunction, np.ones(self._N) )
		alphai = self._inducedAngles(gamma)		
		cl = gamma / self._AREA
		return cl, alphai, self._Y 
