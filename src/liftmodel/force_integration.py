from math import pi as PI
import numpy as np
from scipy.integrate import trapezoid

# module
from src.liftmodel.lifting_line import LiftingLine

class ForceIntegration:
	def __init__(self, N, chord_ratio, twist, advance_ratio=0):
		self._advance_ratio = advance_ratio
		self._lifting_line = LiftingLine(N, chord_ratio, self._geomAngle(twist), self._velRatio)
		
	# --- Private
		
	def _velAngle(self, x):
		return np.arctan( self._advance_ratio * np.abs(x) )	
				
	def _velRatio(self, x):
		return np.sqrt( 1 + (self._advance_ratio * x)**2 ) 
		
	def _geomAngle(self, twist):
		return lambda x: twist(x) - self._velAngle(x)
		
	# --- Public
		
	def setTwist(self, twist, advance_ratio):
		self._advance_ratio = advance_ratio
		self._lifting_line.setGeomAngle( self._geomAngle(twist) )
		self._lifting_line.setVelRatio( self._velRatio )
		
	def integrateCoefficients(self):		
			# solve lift
		[gamma, alphai, y] = self._lifting_line.solveCirculation()
			# airflow angle
		theta = self._velAngle(y) + alphai	
		cos_theta = np.cos(theta)
		sin_theta = np.sin(theta)	
			# quadrature samples
		vel_ratio = self._velRatio(y)
		gamma_vel = gamma * vel_ratio
		gamma_x = gamma_vel * -sin_theta
		gamma_y = gamma_vel * cos_theta
		gamma_z = gamma_y * np.abs(y)	
			# force coefficients
		cy = trapezoid(gamma_y, y)          # force = 0.5 * p * Vinf^2 * Area * coeff
		cx = trapezoid(gamma_x, y)
		ct = trapezoid(gamma_z, y) * 0.5    # torque = 0.5 * p * Vinf^2 * Area * span * coeff
		return cy, cx, ct, gamma, alphai, y
