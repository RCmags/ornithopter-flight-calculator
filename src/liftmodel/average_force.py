from scipy.optimize import fsolve
from math import pi as PI
import numpy as np

# module
from src.liftmodel.force_integration import ForceIntegration

class AverageForce:
	def __init__(self, n, chord_ratio, twist_glide, amplitude, dihedral=0):
		self._force = ForceIntegration(n, chord_ratio, twist_glide)
		self.setAverageY(amplitude, dihedral)
			# store glide data
		data = self._force.integrateCoefficients()[0:5]
		self._glide_coeff = data[0:3]
		self._gamma_glide = data[3]
		self._alpha_glide = self._effectiveAngle( data[4] )		
		
	# --- Private
		
	def _effectiveAngle(self, alphai):
		alpha = self._getGeomAngles() - alphai 
		return alpha
		
	def _getGeomAngles(self):
		return self._force.lifting_line.getGeomAngle()

	# --- Public
	
	def setAverageY(self, amplitude, dihedral):
		self._AVERAGE_Y = np.cos(dihedral) + 0.5*np.cos(0.5*amplitude) - 0.5	
				
	def solveCoefficients(self, twist_up, twist_down, advance_ratio):
			# upstroke
		self._force.setTwist(twist_up, advance_ratio)
		[cy_up, cx_up, ct_up, gamma_up, alphai_up] = self._force.integrateCoefficients()[0:5]
		alpha_up = self._effectiveAngle(alphai_up)
			# downstroke
		self._force.setTwist(twist_down, -advance_ratio)
		[cy_down, cx_down, ct_down, gamma_down, alphai_down] = self._force.integrateCoefficients()[0:5]		
		alpha_down = self._effectiveAngle(alphai_down)
			# averages
		cy_av = 0.5*(cy_up + cy_down) * self._AVERAGE_Y
		cx_av = 0.5*(cx_up + cx_down)	
		ct_max = ct_down if abs(ct_down) > abs(ct_up) else ct_up
			# Store values
		self._average_coeff = [cy_av, cx_av, ct_max]
		self._alpha_up   = alpha_up
		self._alpha_down = alpha_down
		self._gamma_up   = gamma_up	
		self._gamma_down = gamma_down
		
	def solveAdvanceRatio(self, cy, cx, twist_up, twist_down):
		def rootFunc(x):
			[adv_ratio, alpha_root] = x
			fup = lambda x: twist_up(x) + alpha_root
			fdown = lambda x: twist_down(x) + alpha_root                           # pre-calculate twist values to solver faster
				# find averages		
			self.solveCoefficients(fup, fdown, adv_ratio)
			cy_av, cx_av = self._average_coeff[0:2]
			return [cy_av - cy, cx_av - cx]	
		# adv_ratio, alpha_root
		return fsolve( rootFunc, [1, 0] )
		
	def solveAdvanceRatio_LD(self, ld_glide, twist_up, twist_down):
		cy_glide, cx_glide = self._glide_coeff[0:2]
		cx_parasitic = cy_glide / ld_glide + cx_glide                               # counter parasitic drag from body
		return self.solveAdvanceRatio(cy_glide, cx_parasitic, twist_up, twist_down)

	# - Getter functions
		
	def getYaxis(self):
		return self._force.lifting_line.getSpanAxis()
				
	def getChord(self):
		return self._force.lifting_line.getChord()
				
	def getCirculations(self):
		return self._gamma_glide, self._gamma_up, self._gamma_down

	def getAngles(self):
		return self._alpha_glide, self._alpha_up, self._alpha_down
				
	def getGlideCoefficients(self):
		return self._glide_coeff

	def getAverageCoefficients(self):
		return self._average_coeff
