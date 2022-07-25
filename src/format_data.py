from math import pi as PI
import numpy as np
import matplotlib.pyplot as plt

# module
from src.liftmodel.power_balance import PowerBalance
from src.spring_diagram import plotSpring

#---- Constants -----

N_NAME = 22
N_NUM  = 15
N_UNIT = 5
NCHAR  = N_NAME + N_NUM + N_UNIT 

FORMAT_SUBTITLE = ('{:<' + '%d' + '}') % NCHAR
FORMAT_TITLE    = ('{:^' + '%d' + '}') % NCHAR
FORMAT_ROW      = ('{:' + '%d' + '}{:' + '%d' + '}{:' + '%d' + '}') % (N_NAME, N_NUM, N_UNIT)

UNIT_OUTPUT = {'FLIGHT': ['W', 'W', 'W', 'W', '%', '%', 'm', 'm', 'hz'\
                         ,'m/s', 'N.m', 'N.m', 'N', 'N', 'deg']       \
              ,'MOTOR' : ['%', '%', '%', 'W', 'W', 'A', '%']          \
              ,'CLIMB' : ['m/s', 'deg']                               \
              ,'SPRING': ['N/m', 'N', 'm', 'deg'] }

UNIT_INPUT = {'FLIGHT': ['%', 'deg', 'deg', 'm/s^2', 'kg', 'kg', 'kg/m^3', 'm^2'] \
             ,'MOTOR' : ['V', 'A', 'A', 'ohm', 'rpm/V', '%']                      \
             ,'SPRING': ['m', 'm', 'm'] }

#---- functions ------

def labelInputs(n_step,   chord_ratio, twist_glide, twist_up, twist_down                               \
               ,ld_ratio, amplitude,   dihedral,    gravity,  mass_total, mass_wing, air_density, area \
               ,voltage,  i_stall,     i_noload,    motor_rm, motor_kv,   throttle                     \
               ,xoffset,  yoffset,     rspring):
    # create dict               
	input = {}	
	input['WING'] = {                        \
         	 'n_step'          : n_step      \
	    	,'chord_ratio'     : chord_ratio \
	    	,'twist_glide'     : twist_glide \
	    	,'twist_upstroke'  : twist_up    \
	    	,'twist_downstroke': twist_down  }
	
	input['FLIGHT'] = {
	     	 'lift_drag_ratio' : ld_ratio    \
	    	,'amplitude'       : amplitude   \
	    	,'dihedral'        : dihedral    \
	    	,'gravity'         : gravity     \
	    	,'mass_total'      : mass_total  \
	    	,'mass_wing'       : mass_wing   \
	    	,'air_density'     : air_density \
	    	,'area'            : area        }
	    	
	input['MOTOR'] = {                       \
	     	 'voltage'         : voltage     \
	    	,'current_stall'   : i_stall     \
	    	,'current_noload'  : i_noload    \
	    	,'motor_resistance': motor_rm	 \
	    	,'motor_kv'        : motor_kv    \
	    	,'throttle'        : throttle    }
	
	input['SPRING'] = {                      \
	     	 'xoffset'         : xoffset     \
	    	,'yoffset'         : yoffset     \
	    	,'radius_spring'   : rspring     }
	return input

def aircraftCalculations(wing, flight, motor, spring):
	model = PowerBalance(              \
	          wing['n_step']           \
            , wing['chord_ratio']      \
            , wing['twist_glide']      \
            , wing['twist_upstroke']   \
            , wing['twist_downstroke'] \
            , flight['lift_drag_ratio']\
            , flight['amplitude']      \
            , flight['dihedral']       )
	
	data = model.solveDesign(         \
	          spring['xoffset']       \
	        , spring['yoffset']       \
	        , spring['radius_spring'] \
	        , flight['gravity']       \
	        , flight['mass_total']    \
	        , flight['mass_wing']     \
	        , flight['air_density']   \
	        , flight['area']          \
	        , motor['voltage']        \
	        , motor['current_stall']  \
	        , motor['motor_kv']       \
	        , motor['throttle']       \
	        , motor['current_noload'] \
	        , motor['motor_resistance'] )
	
	output = {'FLIGHT'      : data[0] \
	         ,'MOTOR'       : data[1] \
	         ,'CLIMB'       : data[2] \
	         ,'SPRING'      : data[3] \
	         ,'DISTRIBUTION': data[4] }
	return output

def formatNumber(val):
	if abs(val) < 1e-2:
		return format(val, ".3e")
	else:
		return format(val, ".3f")

def printTable(sub_title, data, units):
	# format
	key = list( data.keys() )
	val = list( data.values() )
	val = np.vectorize(formatNumber)(val)
	data = np.transpose( [key, val, units] )
	# display
	print(FORMAT_SUBTITLE.format(sub_title))
	print('-' * NCHAR)
	for row in data:
		print(FORMAT_ROW.format(*row))

def printData(title, data, unit):
	# header
	print('=' * NCHAR)
	print(FORMAT_TITLE.format(title))
	print('=' * NCHAR)
	# tables
	for key in unit.keys():
		printTable( key, data[key], unit[key] )
		print()

def plotLift(x, aglide, aup, adown, glide, up, down):
	fig, axis1 = plt.subplots()
	# lift
	axis1.plot(x, up   , color='blue' , label='lift upstroke')
	axis1.plot(x, down , color='red'  , label='lift downstroke')
	axis1.plot(x, glide, color='black', label='lift glide')	
	# alpha
	axis2 = axis1.twinx() 
	axis2.plot(x, np.degrees(aup)   , color='blue' , alpha=0.6, lw=2 , linestyle='dotted')
	axis2.plot(x, np.degrees(adown) , color='red'  , alpha=0.6, lw=2 , linestyle='dotted')
	axis2.plot(x, np.degrees(aglide), color='black', alpha=0.6, lw=2 , linestyle='dotted')
	# labels
	fig.set_tight_layout(True)
		# axis 1
	axis1.legend(loc='lower center')
	axis1.set_title("lift distribution")
	axis1.set_ylabel("lift per span (N/m)")	
	axis1.set_xlabel("spanwise position (m)")
	axis1.spines['top'].set_visible(False)
	axis1.ticklabel_format(useMathText=True, scilimits=(-2,2))
	axis1.minorticks_on()
		# axis 2
	axis2.set_ylabel("angle of attack (deg)")
	axis2.spines['top'].set_visible(False)
	axis2.minorticks_on()
	
def plotWingGeometry(x, chord, tglide, tup, tdown): 
	fig, axis = plt.subplots(2,1, sharex=True)
	# chord
	y1 = chord * 1/4
	y2 = chord * -3/4
	axis[1].fill_between(x, y1, y2, color='gray', alpha=0.1)
	axis[1].plot(x, y1, color='black', lw=1, label='LE')
	axis[1].plot(x, y2, color='gray' , lw=1, label='TE')	
		# center lines	
	imid = round(len(x) * 0.5) 
	axis[1].vlines(0, y2[imid], y1[imid], color='gray', linestyle='dashed', lw=0.5)
	axis[1].hlines(0, x[0]    , x[-1]   , color='gray', linestyle='dashed', lw=0.5)
	# twist
	axis[0].plot(x, np.degrees(tglide), color='black', label='glide')
	axis[0].plot(x, np.degrees(tup)   , color='blue' , label='upstroke')
	axis[0].plot(x, np.degrees(tdown) , color='red'  , label='downstroke')
	# labels
	fig.set_tight_layout(True)
		# axis 1
	cmax = max(chord)
	axis[1].set_ylim( -cmax, 0.5*cmax )
	axis[1].set_aspect('equal')	
	axis[1].set_ylabel("chord (m)")	
	axis[1].set_xlabel("spanwise position (m)")
	axis[1].spines['top'].set_visible(False)
	axis[1].spines['right'].set_visible(False)
	axis[1].ticklabel_format(useMathText=True, scilimits=(-2,2), axis='y')
	axis[1].minorticks_on()
		# axis 2
	axis[0].legend(loc='lower center')	
	axis[0].set_ylabel("twist angle (deg)")
	axis[0].spines['top'].set_visible(False)
	axis[0].spines['right'].set_visible(False)
	axis[0].ticklabel_format(useMathText=True, scilimits=(-2,2))	
	axis[0].set_title("wing geometry")
	axis[0].minorticks_on()

def inputToRadian(input):
	input['FLIGHT']['dihedral']  = np.radians( input['FLIGHT']['dihedral'] )
	input['FLIGHT']['amplitude'] = np.radians( input['FLIGHT']['amplitude'] )
	# convert functions
	fg = input['WING']['twist_glide']
	fu = input['WING']['twist_upstroke']
	fd = input['WING']['twist_downstroke']
	input['WING']['twist_glide']      = lambda x: np.radians( fg(x) ) 
	input['WING']['twist_upstroke']   = lambda x: np.radians( fu(x) )
	input['WING']['twist_downstroke'] = lambda x: np.radians( fd(x) )
	return input
	
def outputToDegree(input, output):
	input['FLIGHT']['dihedral']  = np.degrees( input['FLIGHT']['dihedral'] )
	input['FLIGHT']['amplitude'] = np.degrees( input['FLIGHT']['amplitude'] )
	output['FLIGHT']['pitch_trim_angle'] = np.degrees( output['FLIGHT']['pitch_trim_angle'] )
	output['CLIMB']['climb_angle']       = np.degrees( output['CLIMB']['climb_angle'] )
	output['SPRING']['spring_angle']     = np.degrees( output['SPRING']['spring_angle'] )
	return input, output

def formatData(input):          # convert inputs and outputs to degrees
	input = labelInputs(*input)
	input = inputToRadian(input)
	output = aircraftCalculations( *input.values() )
	dist = output['DISTRIBUTION']
	# plots	 
	plotWingGeometry( dist['y_axis']   \
          	, dist['chord']            \
          	, dist['twist_glide']      \
          	, dist['twist_upstroke']   \
          	, dist['twist_downstroke'] )
	
	plotLift( dist['y_axis']          \
          	, dist['alpha_glide']     \
          	, dist['alpha_upstroke']  \
          	, dist['alpha_downstroke']\
          	, dist['lift_glide']      \
          	, dist['lift_upstroke']   \
          	, dist['lift_downstroke'] )
	plotSpring( *input['SPRING'].values()           \
              , input['FLIGHT']['dihedral']         \
              , input['FLIGHT']['amplitude']        \
              , output['SPRING']['spring_constant'] )
	# tables
	input, output = outputToDegree(input, output)
	printData("INPUTS", input, UNIT_INPUT)
	printData("OUTPUTS", output, UNIT_OUTPUT)
	plt.show()
