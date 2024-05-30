import numpy as np
import matplotlib.pyplot as plt
# from sympy import symbols, Eq, solve, sqrt

class Ion_Engine:
    """
    Description
    -----------
    - Parametric model of a monoenergetic, monopropellant ion engine 
    - Performance model accounts for ion beam properties, efficiencies and correction factors
    - System decomposition into Propulsion Module, Electrical Module, and Interface Module.
      Interface only applies for SEP configuration.
    - System mass is assessed for the Ion Engine itself and for the Propulsion System in general,
      which is the sum of the 3 modules in the previous point   
      
    Sources
    -------
    [1] Fundamentals of Electric Propulsion: Ion and Hall Thrusters
        D. M. Goebel, I. Katz
        2008
    [2] A Detailed Model of Ion Propulsion Systems
        J. R. Brophy, G. Aston
        1989
    [3] An improved PR equation of state for CO2-containing gas compressibility factor calculation
        L. Zhuoran, J. Wenlong, L. Changjun

    Input Parameters
    ----------------
    thrust (N): float   
        Ion Engine produced thrust
    specific_impulse (s): float 
        Ion Engine specific impulse
    propellant (.): str
        Type of propellant, among the following possibilities of the propellant dictionary
        propellant = {
                      'Xe',
                      'I',
                      }
    design (.): str
        Design of the Ion Thruster, among the following possibilities of the design dictionary
        design = {
                 'ring_cusp',
                 'div_field',
                 }
    power_source (.): str 
        Type of Power source, among the following possibilities of the power_source dictionary
        power_source = {
                 'solar',
                 'nuclear',
                 }
    beam_div_angle (°): float (0 to 90)
        Divergence of the ion beam
    beam_diameter (cm): float
        Diameter of the ion beam
        Default to 28.5 [2]
    beam_ratio (.): float (0 to infinity)
        Ratio of doubly charged current beam over singly charged current beam (0 = no double charge, infinity (put a large number) = no single charge)
    beam_voltage (V): float 
        ion Beam Voltage 
    engine_number (.) = int (min 1)
        Number of Engines in the system
    mass_propellant (kg): float
        Theoretical mass of propellant needed to perform the mission, but the code accounts for residuals
    mass_battery (kg): float
        Mass of the battery system derived from mission analysis. Applies only for Solar Electric Propulsion (SEP)
    mass_solar_array (kg): float
        Mass of the solar array derived from mission analysis. Applies only for Solar Electric Propulsion (SEP)
    """
    def __init__(
                self, 
                thrust,
                specific_impulse,
                propellant,
                design,
                power_source,
                beam_div_angle,
                beam_diameter,
                beam_ratio,
                beam_voltage, 
                engine_number, 
                mass_propellant,
                mass_battery_system,
                mass_solar_array          
                ):
        
        self.thrust = thrust
        self.specific_impulse = specific_impulse
        self.propellant = propellant
        self.design = design
        self.power_source = power_source
        self.beam_div_angle = beam_div_angle 
        self.beam_diameter = beam_diameter
        self.beam_ratio = beam_ratio 
        self.beam_voltage = beam_voltage
        self.engine_number = engine_number
        self.mass_propellant = mass_propellant
        self.mass_battery_system = mass_battery_system
        self.mass_solar_array = mass_solar_array
            
            
        """ Constants """
        
        
        earth_gravity = 9.807 # m/s2, Earth's gravity at sea-level
        elem_charge = 1 # eV, Elementary charge (= electron charge)
        eV_to_J = 1.602176487e-19 # J, Conversion from eV to J
        Ion_mass = { # AMU = g/mol 
            'Xe': 131.29,
            'I': 126.90,
            }
        critical_pressure = { # (MPa) 
            'Xe': 5.841,
            'I': 6.279,
            } 
        critical_temperature = { # (K) 
            'Xe': 289.74,
            'I': 819.0,
            } 
        acentric_factor = { # (.) 
            'Xe': 0.008,
            'I': 0.229,
            }               
        AMU_to_kg = 1.6605e-27 # Conversion from AMU to kg
        gas_constant = 8.314 # J/(mol.K)
        
        
        """ Thrust Correction Factors """
        
        
            # Beam Divergence
        
        div_corr = np.cos(np.deg2rad(beam_div_angle)) # (.), takes into account beam diverging exhaust, Eq. (2.3-10) [1] 
        
            # Doubly Charged Ions
        
        eff_thrust_double_charge = (1 + (1/np.sqrt(2))*beam_ratio)/(1 + beam_ratio) # (.), takes into account beams of singly and doubly charged ions, Eq. (2.3-14) [1] 
    
            # Thrust Total Correction
            
        thrust_corr = div_corr*eff_thrust_double_charge # (.), Eq. (2.3-15) [1]
        
        
        """ Engine Operating Point """
        
        
            # Ion Beam Current & Voltage
        
        self.second_member = (thrust/thrust_corr)*np.sqrt(elem_charge*eV_to_J/(2*Ion_mass[propellant]*AMU_to_kg)) # Eq. (2.3-16)
        self.beam_current = self.second_member/np.sqrt(beam_voltage) # A, Eq. (2.3-16) [1]
        
            # Mass Flow Rates
                
        mass_flow_ion = self.beam_current*Ion_mass[propellant]*AMU_to_kg/(elem_charge*eV_to_J) # kg/s, Eq. (2.3-7) [1]
        self.mass_flow_prop = thrust/(earth_gravity*specific_impulse) # kg/s, Eq. (2.4-1) [1]
        
            # Mass Efficiency
                
        eff_mass_double_charge = (1 + (1/2)*beam_ratio)/(1 + beam_ratio) # (.), TO BE CHECKED takes into account beams of singly and doubly charged ions, Eq. (2.4-7) [1]
        efficiency_mass = eff_mass_double_charge*(mass_flow_ion/self.mass_flow_prop)
    
        
        self.beam_voltage_min = ((self.second_member*Ion_mass[propellant]*AMU_to_kg*specific_impulse*earth_gravity)/(elem_charge*eV_to_J*thrust))**2 # V, mass_flow_ion must be less than mass_flow_prop 
        
        if efficiency_mass > eff_mass_double_charge:
            raise ValueError(f'Not realistic as is, the voltage must the greater than: {self.beam_voltage_min} V')
                          
        
        """ Thruster Efficiency """
        
        
        # Total Efficiency
        
        efficiency_electrical_guess = 0.85
        eff_thruster_total_guess = (thrust_corr**2)*efficiency_electrical_guess*efficiency_mass # (.), Eq. (2.5-7) [1], to be updated in the PPU loop
        
    
        """ System Mass Decomposition, Cf [2] for equations that are not numbered (nor are they in the paper) """
        
        
        # --- ELECTRICAL POWER SYSTEM (Mass_EPS) ---
        
        
            # Power Processing Unit
        
        flight_packaging_factor = 2.0 # between 1.7 and 2.0
        
        PPU_efficiencies = { # (.) 
            'filter': 0.99,
            'beam': 0.95,
            'discharge': 0.88,
            'inverter': 0.97,
            'LLS': 0.70,
            } 
        
        mass_inverter = 0.75 # kg, mass of DC/AC inverter assumed constant 
        mass_LLS = 1.50 # kg, mass of Low Level Supplies assumed constant 
        mass_CCT = 0.40 # kg, mass of Command/Control/Telemetry assumed constant 
        
        
        """ Power Breakdown """
        
        
            # Jet Power
        
        power_jet = (thrust**2)/(2*self.mass_flow_prop)*10**(-3) # kW, Eq. (2.5-4) [1], kinetic power of the jet 
        
            # Total Power
        
        power_engine_guess = power_jet/eff_thruster_total_guess # kW, Eq. (2.5-3) [1], total input power needed to make the spacecraft move as expected. Equal to power_beam/efficiency_electrical if beam_ratio = 0      
        k = 0
        eps = 0.00001
        
            # Iterative loop to determine the final efficiencies & input power

        power_beam_output = self.beam_current*beam_voltage*10**(-3) # kW, Eq. (2.5-1) [1]
        power_beam_input = power_beam_output/PPU_efficiencies['beam'] # kW
        power_discharge_output = power_engine_guess - power_beam_output - 50*10**(-3)*PPU_efficiencies['inverter']*PPU_efficiencies['LLS'] # kW, ERROR in paper (P_BO instead of P_BI) + inverter input power is assumed constant at 50 W
        power_discharge_input = power_discharge_output/PPU_efficiencies['discharge'] # kW
        power_filter_output = power_beam_input + power_discharge_input + 50*10**(-3) # kW, ERROR in paper (input powers and not outputs)
        power_filter_input = power_filter_output/PPU_efficiencies['filter'] # kW
        efficiency_PPU = power_engine_guess/power_filter_input 
        self.eff_thruster_total = (thrust_corr**2)*efficiency_PPU*efficiency_mass # (.), Eq. (2.5-7) [1]
        
        while abs(eff_thruster_total_guess - self.eff_thruster_total) > eps and k < 1000:
            
            eff_thruster_total_guess = self.eff_thruster_total
            self.power_engine = power_jet/self.eff_thruster_total # kW
            power_discharge_output = self.power_engine - power_beam_output - 50*10**(-3)*PPU_efficiencies['inverter']*PPU_efficiencies['LLS'] # kW
            power_discharge_input = power_discharge_output/PPU_efficiencies['discharge'] # kW
            power_filter_output = power_beam_input + power_discharge_input + 50*10**(-3) # kW
            power_filter_input = power_filter_output/PPU_efficiencies['filter'] # kW
            efficiency_PPU = self.power_engine/power_filter_input # (.)
            self.eff_thruster_total = (thrust_corr**2)*efficiency_PPU*efficiency_mass # (.)
            k = k + 1
            #print(k)
                   
            # PPU Heat dissipation
            
            power_dissipation_PPU = power_filter_input - self.power_engine # kW, from [2] but equvalent to Eq. (2.6-1) [1], power that needs to be radiated away
            self.power_heat_engine = power_dissipation_PPU
        
            # Thrust to Power Ratio
            
        self.T_to_P_ratio = thrust*10**3/(self.power_engine) # mN/kW, Eq. (2.5-9) [1]
        
        mass_filter = 1.03 + 0.502*power_filter_output # kg
        mass_beam_supply = 2.03 + 1.43*power_beam_output # kg
        mass_discharge_supply = 2.03 + 1.43*power_discharge_output # kg
        
        mass_PPU = engine_number*flight_packaging_factor*(mass_inverter + mass_LLS + mass_CCT + mass_filter + mass_beam_supply + mass_discharge_supply) # kg
        
            # Power Processing Thermal Control (Assumes maximum PPU base plate temperature of 50°C)
        
        mass_pipes = 12.5*power_dissipation_PPU # kg, stainless steel/alcohol heat pipes
        mass_radiator = 8.0*power_dissipation_PPU # kg, aluminium radiators
        
        mass_PPTC = engine_number*(mass_pipes + mass_radiator) # kg
        
            # Power Subsystem Thermal Control
            
        power_dissipation_nonPPU = (0.4/6)*self.engine_number # kW, [2] 
        mass_sub_thermal = 31.8*power_dissipation_nonPPU # kg
        
            # Total Heat Generated
        
        self.power_dissipation = self.engine_number*power_dissipation_PPU + power_dissipation_nonPPU # kW, approximation given that the design point in [2] is for 6 engines        
            
            # High Voltage Power Distribution 
        
        mass_HVPD = (0.5 + 0.5)*engine_number # kg, one switch/fuse per PPU = per engine, + accounts for wiring
        
            # DC/DC Converter
        
        converter_redudancy = 1
        mass_converter = (1 + converter_redudancy)*flight_packaging_factor*(2.03 + 1.43*power_discharge_output) # kg, same scaling relationship as PPU discharge supply (error in paper)
        
            # Low Voltage Power Distribution 
        
        mass_distrib_assembly = 1.48*engine_number # kg
        mass_pyro_switch = 0.38*engine_number # kg   
        mass_power_switch = 2.2 # kg, accounts for redudancy
        
        mass_LVPD = mass_distrib_assembly + mass_pyro_switch + mass_power_switch # kg
            
            # Data Handling System
        
        mass_DHS = 2*engine_number # kg
     
            # Battery System, has to be added for solar configuration. Battery sizing comes from mission design (eclipse time)
        
            # Power Subsystem Structure
        
        if power_source == 'nuclear':
            mass_PSS = 0.08*(mass_PPU + mass_PPTC + mass_sub_thermal + mass_converter + mass_HVPD + mass_LVPD + mass_DHS) # kg
        elif power_source == 'solar':
            mass_PSS = 0.08*(mass_PPU + mass_PPTC + mass_sub_thermal + mass_converter + mass_HVPD + mass_LVPD + mass_DHS + mass_battery_system) # kg

            # Total Electrical Power System Mass
            
        if power_source == 'nuclear':
            self.Mass_EPS = mass_PPU + mass_PPTC + mass_sub_thermal + mass_HVPD + mass_converter + mass_LVPD + mass_DHS + mass_PSS # kg
        elif power_source == 'solar':
            self.Mass_EPS = mass_PPU + mass_PPTC + mass_sub_thermal + mass_HVPD + mass_converter + mass_LVPD + mass_DHS + mass_PSS + mass_battery_system # kg
            
        
        # --- PROPULSION MODULE (Mass_PM) ---
        
        
            # Engine
            
        if design == 'ring_cusp': # Eq for M_e, p.4 [2] 
            self.mass_engine = -0.382 + 0.424*beam_diameter + 0.00297*beam_diameter**2 # kg
        elif design == 'div_field':
            self.mass_engine = 0.346 + 0.199*beam_diameter + 0.00288*beam_diameter**2 # kg
    
            # Gimbal system
        
        mass_gimbal = 0.3*self.mass_engine # kg, Eq for M_gimbal, p.4 [2], assuming a conventional two-axis gimble ring configuration, and the fraction remains constant for all engines
        
            # Tank (single tank) 
          
        # It requires more info on how to derive compressibility factor for Iodine 'I', NOT TRIVIAL
        
        # Based on: An improved PR equation of state for CO2-containing gas compressibility factor calculation, Li. Z, Jia. W, Li. C, 2016
        
        # reduced_temperature = tank_temperature/critical_temperature[propellant]   
        # alpha_coeff = (1 + (1 - np.sqrt(reduced_temperature))*(0.37464 + 1.54226*acentric_factor[propellant] - 0.26992*acentric_factor[propellant]**2))**2
        # a = 0.45724*alpha_coeff*((gas_constant*critical_temperature[propellant])**2)/(critical_pressure[propellant]*10**6)
        # b = 0.07780*gas_constant*critical_temperature[propellant]/(critical_pressure[propellant]*10**6)
        
        # A = a*tank_pressure*10**6/((gas_constant*tank_temperature)**2)
        # B = b*tank_pressure*10**6/(gas_constant*tank_temperature)
        # w = 1
        # x = -(1 - B)
        # y = A - 2*B - 3*B**2
        # z = -(A*B - B**2 - B**3)        
        # coefficients = [w, x, y, z]
        # roots = np.roots(coefficients)
        # print(roots)
        # compressibility_factor =  [root for root in roots if 0 <= root <= 1]
        # print(compressibility_factor)
        # propellant_density = tank_pressure*(10**6)*Ion_mass[propellant]*10**(-3)/(compressibility_factor*gas_constant*tank_temperature) # kg/m3, Van der Waals law
        # print(propellant_density)
    
        mass_usable = 1.04*mass_propellant # kg
        mass_tank = 0.137*mass_usable # kg
        mass_blankets = 0.01*mass_usable # kg
        mass_residual = 0.017*mass_usable # kg
        mass_total_propellant = mass_usable + mass_residual # kg
        
        mass_tank = mass_tank + mass_blankets # kg
               
        # tank_radius = (3*mass_total_propellant/(4*np.pi*propellant_density))**(1/3)
        
            # Propellant Distribution Systems
        
        mass_PDS = 4.45 + 2.26*engine_number # kg
        
            # Propulsion Module Cabling
        
        mass_GWB = 0.27*engine_number # kg, mass of gimbal system wiring
        mass_PDSC = 0.17*engine_number # kg, mass of propellant system wiring 
        mass_EW = 0.16*engine_number*self.power_engine # kg, mass of PPU to engine wiring
        
        mass_PMC = mass_GWB + mass_PDSC + mass_EW # kg
        
            # Propulsion Module Structure
        
        mass_PMS =  0.04*(mass_tank + mass_total_propellant + mass_PDS + mass_GWB + mass_EW + mass_PDSC + engine_number*(self.mass_engine + mass_gimbal)) # kg, TO BE CHECKED, ref [2] says 0.4 but should be 0.04
        
        # Total Propulsion Module Mass
        
        self.Mass_PM = engine_number*(self.mass_engine + mass_gimbal) + mass_tank + mass_PDS + mass_PMC + mass_PMS # kg
        
        
        # --- POWER INTERFACE SYSTEM (Mass_PIS) ---
        
        
        mass_roll_ring = engine_number*0.83*self.power_engine # kg
        mass_articulation = 0.015*mass_solar_array # kg
        mass_PMM = 3 # kg, Power Monitor Unit independent from solar array 
        mass_components = mass_roll_ring + mass_articulation + mass_PMM # kg
        mass_PIS_structure = 0.02*mass_solar_array + 0.08*(mass_components) # kg 
        
        # Total Power Interface System Mass
        
        if power_source == 'nuclear':
            self.Mass_PIS = 0
        elif power_source == 'solar':
            self.Mass_PIS = mass_components + mass_PIS_structure # kg

        
        # --- TOTAL MASS --- 
            
            
        self.mass_dry = self.Mass_PM + self.Mass_EPS + self.Mass_PIS # kg
        self.mass_wet = self.mass_dry + mass_total_propellant # kg


if __name__ == "__main__":

    Ion_Engine_instance = Ion_Engine(
        thrust = 207*10**(-3), # N
        specific_impulse = 3500, # s
        propellant = 'Xe', # 'Xe' or 'I', both FINE for engine performance, but WARNING! requires further investigation for tank calculations with 'I'
        design = 'ring_cusp', # 'ring_cusp' or 'div_field'
        power_source = 'solar', # 'solar' or 'nuclear'
        beam_div_angle = 10, # 0° to 90° 
        beam_diameter = 28.5, # cm
        beam_ratio = 0.1, # 0 to infinity
        beam_voltage = 1205, # V
        engine_number = 6,
        mass_propellant = 446.5, # kg, theoretical mass needed, and code accounts for residuals
        mass_battery_system = 48, # kg, applicable if power_source = 'solar'. Can put whatever value otherwise
        mass_solar_array = 400 # kg, applicable if power_source = 'solar'. Can put whatever value otherwise
        )    

    engine_number = Ion_Engine_instance.engine_number
    engine_mass = Ion_Engine_instance.mass_engine
    beam_voltage = Ion_Engine_instance.beam_voltage
    beam_current = Ion_Engine_instance.beam_current
    second_member = Ion_Engine_instance.second_member
    thrust = Ion_Engine_instance.thrust
    specific_impulse = Ion_Engine_instance.specific_impulse
    eff_thruster_total = Ion_Engine_instance.eff_thruster_total
    mass_flow_prop = Ion_Engine_instance.mass_flow_prop
    T_to_P_ratio = Ion_Engine_instance.T_to_P_ratio
    power_engine = Ion_Engine_instance.power_engine
    power_dissipation = Ion_Engine_instance.power_dissipation
    power_heat_engine = Ion_Engine_instance.power_heat_engine
    mass_PM = Ion_Engine_instance.Mass_PM
    mass_EPS = Ion_Engine_instance.Mass_EPS
    mass_PIS = Ion_Engine_instance.Mass_PIS
    mass_dry = Ion_Engine_instance.mass_dry
    mass_wet = Ion_Engine_instance.mass_wet
    
    print(f'\n Ion Engine: \n')
    
    print(f'Thrust (mN) = {thrust*10**3}')
    print(f'Specific Impulse (s) = {specific_impulse}')
    print(f'Beam Voltage (V) = {beam_voltage}')
    print(f'Beam Current (A) = {beam_current}')
    print(f'Total Thruster Efficiency (.) = {eff_thruster_total}') 
    print(f'Engine Input Power (PPU Output) (kW) = {power_engine}')
    print(f'Thrust to Power ratio (mN/kW) = {T_to_P_ratio}')
    print(f'Engine Heat Power (to be dissipated) (kW) = {power_heat_engine}')
    print(f'Propellant Mass Flow Rate (g/s) = {mass_flow_prop*10**3}')
    print(f'Engine Mass (kg) = {engine_mass}')
    
    print(f'\n Propulsion Module: \n')
    
    print(f'Number of Engine(s) = {engine_number}')
    print(f'Total Thrust (mN) = {engine_number*thrust*10**3}')
    print(f'Total Engines Input Power (kW) = {engine_number*power_engine}')
    print(f'Total Heat Power (to be dissipated) (kW) = {power_dissipation}')
    print(f'PM Mass (kg) = {mass_PM}')
    print(f'EPS Mass (kg) = {mass_EPS}')
    print(f'PIS Mass (kg) = {mass_PIS}')
    print(f'Dry Mass (kg) = {mass_dry}')
    print(f'Wet Mass (kg) = {mass_wet}')
     
   
    # # Derivation of the curve I = f(V) to find the admissible values of beam current and voltage to provide the required thrust and specific impulse
        
    # beam_voltage_plot = np.linspace(0, max(beam_current,beam_voltage), 100)    
    # beam_current_plot = (second_member/np.sqrt(beam_voltage_plot)) 

    # plt.plot(beam_voltage_plot, beam_current_plot, color = 'red')
    # plt.xlabel('Beam Voltage (V)')
    # plt.ylabel('Beam Current (I)')
    # plt.title(f"Beam Current as a function of Beam Voltage \n for Thrust = {thrust*10**3:.2f} mN and Specific Impulse = {specific_impulse} s")
    # plt.grid(True)
    # plt.show()