from math import pi as PI
import numpy as np
from scipy.integrate import quad

# module
from src.liftmodel.average_force import AverageForce

class PowerBalance:
	def __init__(self, n, chord_ratio, twist_glide, twist_up, twist_down, ld_ratio, amplitude, dihedral=0):
		# lift model
		self.liftmodel = AverageForce(n, chord_ratio, twist_glide, amplitude, dihedral)		
		self._ADV_RATIO, self._PITCH_TRIM = self.liftmodel.solveAdvanceRatio_LD(ld_ratio, twist_up, twist_down)
		self._LD_RATIO  = ld_ratio
		self._AMPLITUDE = amplitude 
		self._DIHEDRAL  = dihedral
		# inertia
		area_moment_0 = quad(chord_ratio, 0, 1)[0]
		area_moment_2 = quad( lambda x: chord_ratio(x)*(x**2), 0, 1)[0]
		self._CONST_SPAN = 0.5 / area_moment_0
		self._CONST_INERTIA_MASS = area_moment_2 / area_moment_0 
				
		# other
		self._CONST_FREQ = np.sqrt(2)/( PI * amplitude )
		self._CHORD_RATIO_0 = chord_ratio(0)
		# wing twist
		yaxis = self.liftmodel.getYaxis()
		self._TWIST = [ np.vectorize(twist_glide)(yaxis)\
		              , np.vectorize(twist_up)(yaxis)   \
		              , np.vectorize(twist_down)(yaxis) ] 		
	# --- Private
		
	def _flightPower(self, gravity, mass_total, mass_wing, air_density, area):
		# coefficients
		cy_av   = self.liftmodel.getAverageCoefficients()[0]
		ct_peak = self.liftmodel.getAverageCoefficients()[2]
		ct_g    = self.liftmodel.getGlideCoefficients()[2]
		# forces
		lift   = mass_total * gravity
		thrust = lift / self._LD_RATIO
		# wing motion	
		vel     = np.sqrt( 2 * lift / ( air_density * area * cy_av ) )
		radius  = np.sqrt( self._CONST_SPAN * (area * 0.5) )
		ang_vel = self._ADV_RATIO * vel / radius
		chord   = self._CHORD_RATIO_0 * radius * 2 
		# inertia loss
		inertia = self._CONST_INERTIA_MASS * mass_wing * radius**2
		freq    = ang_vel * self._CONST_FREQ
		# aerodynamic torques per half span
		torque       = (lift/cy_av) * radius
		torque_glide = torque * ct_g
		torque_max   = torque * ct_peak
		# power draw
		power_min      = thrust * vel
		power_inertial = inertia * freq * ang_vel**2 
		power_aerodyn  = 2 * (torque_max - torque_glide) * ang_vel	# total torque for two wings
		power_mech     = power_inertial + power_aerodyn
		# efficiencies
		eff_aero = power_min / power_aerodyn
		eff_mech = power_min / power_mech
		
		output = {'power_mechanical' : power_mech    \
		         ,'power_aerodynamic': power_aerodyn \
		         ,'power_inertial'   : power_inertial\
		         ,'power_minimum'    : power_min     \
		         ,'eff_aerodynamic'  : eff_aero      \
		         ,'eff_mechanical'   : eff_mech      \
		         ,'wing_radius'      : radius	     \
		         ,'wing_root_chord'  : chord         \
		         ,'frequency'        : freq          \
		         ,'velocity'         : vel 	         \
		         ,'torque_glide'     : torque_glide  \
		         ,'torque_max'       : torque_max    \
		         ,'force_lift'       : lift          \
		         ,'force_thrust'     : thrust        \
		         ,'pitch_trim_angle' : self._PITCH_TRIM}
		return output
		
	def _motorPower(self, power, vin, istall, kv, freq_wing, i_noload=0):
		# max available power
		power_max = vin * istall * 0.25	
		
		# loading factors		
		cp = power / power_max					# power ratio
		if cp > 1:
			cw = 0.5
			print("WARNING: Insufficient power to maintain level flight \n")
		else:
			cw = 0.5 + 0.5*np.sqrt(1 - cp)		# angular velocity
			
		# electrical power
		iload     = istall*(1 - cw)
		resistor  = vin / istall
		power_in  = power + i_noload*vin + resistor*iload**2
		imotor    = power_in / vin
		eff_motor = power / power_in
 		# motor speed 
		freq_max   = vin * kv / 60		
		freq_motor = freq_max * cw
		gear_ratio = freq_motor / freq_wing
		
		output = {'resistance'          : resistor     \
		         ,'coeff_power'         : cp           \
		         ,'coeff_angvel'        : cw           \
		         ,'efficiency_motor'    : eff_motor    \
		         ,'power_electric_input': power_in     \
		         ,'power_max_mechanical': power_max    \
		         ,'current_motor'       : imotor       \
		         ,'gear_ratio'          : gear_ratio}
		return output	

	def _climbPower(self, weight, vel, power_max, power_flight, power_inertia): 
		# power for thrust
			# constant amplitude, higher frequency -> higher inertial loss
		power_inertia_change = power_inertia*( (power_max/power_flight)**(3/2) - 1 )	
		power_excess = power_max - (power_flight + power_inertia_change)
		# climb rate			
		climb_rate  = power_excess / weight			
		climb_angle = np.arcsin(climb_rate / vel)
				
		output = {'climb_rate' : climb_rate \
		         ,'climb_angle': climb_angle} 
		return output

	def _geometry(self, x, y, ro, ang):
		dy = ro*np.sin(ang) + y
		dx = ro*np.cos(ang) - x
		mag = np.sqrt(dx**2 + dy**2)
		return mag, dx, dy

	def _springBias(self, x, y, ro, torque):
		# angles
		beta  = self._DIHEDRAL
		alpha = self._AMPLITUDE * 0.5
		# lengths
		len_mid, dx, dy = self._geometry(x, y, ro, beta)
		len_up          = self._geometry(x, y, ro, beta + alpha)[0]
		len_down        = self._geometry(x, y, ro, beta - alpha)[0]
		# moment arm
		dlen   = len_up - len_down
		dlen_k = len_mid - len_down
		angle  = beta + np.arctan2(dy, dx)
		radius = ro * np.sin(angle)
		# spring values
		force    = torque / radius
		k_spring = force / dlen_k
		
		output = {'spring_constant': k_spring\
		         ,'spring_force'   : force   \
		         ,'length_change'  : dlen    \
		         ,'spring_angle'   : angle   } 
		return output

	def _scaleDistributions(self, density, vel, area, radius):
		# scalars
		span = 2 * radius
		dldy = density * area * vel**2 / span
		# planform
		yaxis = self.liftmodel.getYaxis() * radius
		chord = self.liftmodel.getChord() * span
		# lift per span
		gamma_glide, gamma_up, gamma_down = np.array( self.liftmodel.getCirculations() ) * dldy
		# angles
		alpha_glide, alpha_up, alpha_down = self.liftmodel.getAngles() 
		twist_glide, twist_up, twist_down = self._TWIST
				
		output = {'y_axis'          : yaxis       \
		         ,'chord'           : chord       \
		         ,'lift_glide'      : gamma_glide \
                 ,'lift_upstroke'   : gamma_up    \
                 ,'lift_downstroke' : gamma_down  \
                 ,'alpha_glide'     : alpha_glide \
                 ,'alpha_upstroke'  : alpha_up    \
                 ,'alpha_downstroke': alpha_down  \
                 ,'twist_glide'     : twist_glide \
                 ,'twist_upstroke'  : twist_up    \
                 ,'twist_downstroke': twist_down  }
		return output
	
	# --- Public

	def solveDesign(self, xoffset, yoffset      , lever_spring,                    \
	                      gravity, mass_total   , mass_wing   , air_density, area, \
	                      voltage, current_stall, motor_kv    ,                    \
	                      current_noload=None   ):
	
		flight = self._flightPower( gravity     \
		                          , mass_total  \
		                          , mass_wing   \
		                          , air_density \
		                          , area ) 
		                          
		motor = self._motorPower( flight['power_mechanical'] \
		                        , voltage                    \
		                        , current_stall              \
		                        , motor_kv                   \
		                        , flight['frequency']        \
		                        , current_noload             )

		climb = self._climbPower( flight['force_lift']          \
		                        , flight['velocity']            \
		                        , motor['power_max_mechanical'] \
		                        , flight['power_mechanical']    \
		                        , flight['power_inertial'] )
		
		spring = self._springBias( xoffset                \
		                         , yoffset                \
		                         , lever_spring           \
		                         , flight['torque_glide'] )
		                         
		distribution = self._scaleDistributions( air_density           \
		                                       , flight['velocity']    \
		                                       , area                  \
		                                       , flight['wing_radius'] )
		# output data
		return flight, motor, climb, spring, distribution
		
