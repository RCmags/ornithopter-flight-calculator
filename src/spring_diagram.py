from math import pi as PI
import numpy as np
from matplotlib import pyplot as plt

# constants
MID  = 49			# middle index of angle steps 
NUM  = 1 + 2*MID  	# number of angle steps
FREQ = 20			# number of sine wave peaks
LEN_RATIO = 0.4		# relative length of spring 

# functions

def plotSpringGeometry(axis, x, y, r, dih, amp):
	# coordinates
	amp *= 0.5	
	angle = np.linspace( dih + amp, dih - amp, num=NUM)
	xr = np.cos(angle)*r
	yr = np.sin(angle)*r + y
	
	# magnitude and slope angle
	dx = xr - x
	dy = yr
	mag = np.sqrt(dx**2 + dy**2)
	sangle = np.arctan2(dy, dx)
	
	# sine wave for spring
	mag_m = mag[MID]
	xs = np.linspace(0, 1, num=NUM)
	ys = np.sin(FREQ*PI*xs) * mag_m * 0.05
	ys = [0, *ys, 0]
	
	def rotatePoints(l, angle):
		xr = (l - LEN_RATIO*mag_m)*xs + 0.5*LEN_RATIO*mag_m
		xr = [0, *xr, l]
		# rotate vector
		rot = np.array([ [np.cos(angle), -np.sin(angle)],
      	                 [np.sin(angle),  np.cos(angle)] ])
		xr, yr = np.dot(rot, [xr, ys])
		xr += x
		return xr, yr
	
	# rotate spring coordinates
	xsm, ysm = rotatePoints( mag_m  , sangle[MID] )	# middle
	xst, yst = rotatePoints( mag[0] , sangle[0] )	# top
	xsb, ysb = rotatePoints( mag[-1], sangle[-1] )	# bottom

	# plots	
		# springs
	axis.plot( xsm, ysm, color='grey' ) 
	axis.plot( xst, yst, color='blue', alpha=0.3 )
	axis.plot( xsb, ysb, color='red' , alpha=0.3 )
		# glide lines
	axis.plot(     xr,          yr     , color='grey', alpha=0.6, lw=1, linestyle='dashed')
	axis.plot( [0, xr[0] ], [y, yr[0] ], color='blue', alpha=0.3, lw=1, linestyle='solid' )
	axis.plot( [0, xr[-1]], [y, yr[-1]], color='red' , alpha=0.3, lw=1, linestyle='solid' )
	axis.hlines(  [y]  ,  0  ,  r      , color='grey', alpha=0.6, lw=1, linestyle='dashed')
		# frame
	axis.plot( [0 , 0]      , [0, y]      , color='black', lw=2 )
	axis.plot( [0 , x]      , [0, 0]      , color='black', lw=2 )
	axis.plot( [0 , xr[MID]], [y, yr[MID]], color='black', marker='o' )
	axis.scatter( [x]       , [0]         , color='black' )
	
	# figure settings
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)	
	axis.ticklabel_format(useMathText=True, scilimits=(-2,2))
	axis.set_xlabel('x (m)')
	axis.set_ylabel('y (m)')
	axis.set_aspect('equal')
	return mag, sangle, angle

def plotSpringTorque(axis, r, mag, sangle, angle, kspring):
	# torque about root
	theta  = sangle - angle
	dist   = np.sin(theta)*r 
	force  = ( mag - mag[-1] )*kspring
	moment = force * dist
	
	# plots
		# torque plot
	angle = np.degrees(angle)
	axis.plot( angle, moment, color='black')
		# endpoints
	axis.scatter( [angle[0]]  , [moment[0]]  , color='blue', marker='o', label='expand') 
	axis.scatter( [angle[-1]] , [moment[-1]] , color='red' , marker='o', label='contract')
	axis.scatter( [angle[MID]], [moment[MID]], color='gray', marker='o', label='glide')
		# center line
	axis.hlines( [moment[MID]], angle[0], angle[-1], color='gray', alpha=0.8, linestyle='dashed')
		# regions
	axis.fill_between( angle[MID:-1], moment[MID:-1], moment[MID], color='red' , alpha=0.2 )
	axis.fill_between( angle[0:MID] , moment[0:MID] , moment[MID], color='blue', alpha=0.2 )
	
	# figure settings
	axis.spines['top'].set_visible(False)
	axis.spines['right'].set_visible(False)	
	axis.ticklabel_format(useMathText=True, scilimits=(-2,2))
	axis.set_xlabel('wing angle (deg)')
	axis.set_ylabel('root torque (N.m)')
	axis.legend()

def plotSpring(x, y, r, dih, amp, kspring):
	# set figure
	fig, axis = plt.subplots(1, 2)
	fig.set_tight_layout(True)
	fig.suptitle('spring geometry')
	# plots
	output = plotSpringGeometry(axis[0], x, y, r, dih, amp)	
	plotSpringTorque(axis[1], r, *output, kspring)

