from scipy.integrate import trapezoid
from math import pi as PI
import numpy as np

class LiftingLine:
	def __init__(self, n, chord_ratio, geom_angle):
		self._theta = np.linspace(0, PI, num=n)
		self._nval = np.arange(1, n+1)
		self._x = np.cos(self._theta) 
		# sine map
		fsin = lambda n, theta: np.sin(n*theta)
		self._sint = self._mapFunc(fsin, self._nval, self._theta) 
		# coefficients
		self._setCoeff1()
		self.setChord(chord_ratio)
		self.setGeomAngle(geom_angle)
	
	#----- private
	
	def _mapFunc(self, func, xaxis, yaxis):
		return func(xaxis[:,None], yaxis[None,:]).T

	def _fsin(self, x, y):
		return self._sint[y-1, x-1]

	def _sinRatioA(self, n, theta):
		return n * self._fsin(n, theta) / self._fsin(1, theta)
	
	def _sinRatioB(self, n, itheta):
		return (2/PI) * self._fsin(n, itheta) / self._chord
	
	def _setCoeff1(self):
		ratio = self._mapFunc(self._sinRatioA, self._nval, self._nval[1:-1])
		n_sq  = self._nval**2
		n_odd = np.where( self._nval % 2, n_sq, -n_sq )
		output = np.insert( ratio, 0, n_sq )
		output = np.append( output, n_odd )
		nmax = self._nval[-1] 
		self._k1 = output.reshape(nmax, nmax)
								
	def _setCoeff2(self):
		k2 = self._mapFunc(self._sinRatioB, self._nval, self._nval)
		self._k = self._k1 + k2
		
	#----- public
		
	def setGeomAngle(self, geom_angle):
		self._alpha = np.vectorize(geom_angle)(self._x)

	def setChord(self, chord_ratio):
		chord = np.vectorize(chord_ratio)(self._x)
		self._chord = np.where( chord == 0, 1e-8, chord )
		self._area = -trapezoid(chord, self._x)             # dimensionless area, neg due to cos map
		self._setCoeff2()                                   # update coefficient matrix		
		
	def getGeomAngle(self):
		return self._alpha
	
	def getChord(self):
		return self._chord
		
	def getSpanAxis(self):
		return np.flip(self._x)
	
	def solveCirculation(self):
		An = np.linalg.solve( self._k, self._alpha )		
		alphai = np.sum( An*self._k1, axis=1 )
		gamma = np.sum( An*self._sint, axis=1 )
		cl = gamma * 4/self._area
		return np.flip(cl), np.flip(alphai), self.getSpanAxis() # rotate due to consine map
