import json
import numpy as np

# np.set_printoptions(precision=15)

# Here are some general functions that don't need to be a part of the pump class structure
def fahrenheit_to_kelvin(fahrenheit):
    """This function takes in a temp in fahrenheit and spits out the kelvin equivalent"""
    kelvin = (fahrenheit - 32)*(5/9) + 273.15
    return kelvin


def psi_to_Pa(psi):
    """This function takes in a pressure in psi and spits out the Pa equivalent"""
    Pa = psi*6894.76
    return Pa


def feet_to_m(feet):
    """This function takes in a length in feet and spits out the equivalent in meters"""
    m = feet*0.3048
    return m


def calc_cyl_cross_area(Diam):
        """This function calculates the hydraulic area of the annulus using the following inputs:
        
        Parameters
        ----------
        - Diam: Diameter of a cylinder [m] (float)
        """

        cross_area = (np.pi*Diam**2)/4
        return cross_area


def calc_mdot(rho, v, hyd_area):
         """This function returns a mass flow rate given the following inputs:
            
            Parameters
            ----------
            - rho: density [kg/m^3] (float)
            - v: velocity [m/s] (float)
            - hyd_area: cross-sectional area through which a fluid is flowing [m^2] (float)"""
         mdot = rho*v*hyd_area
         return mdot


def calc_vdot(rho, mdot):
         """This function returns a volumetric flow rate given the following inputs:
            
            Parameters
            ----------
            - rho: density [kg/m^3] (float)
            - mdot: mass flow rate [kg/s] (float)"""
         
         vdot = mdot/rho
         return vdot


class pump:
    """This class contains pump attributes that help select a correct pump for a home brewing kit"""
    def __init__(self, json_file):
        """This loads all the initial values for the iteration processes"""
        self.json_file = json_file
        self.load_json()


    def load_json(self):
        """This function pulls in all the input values from the json"""
        with open(self.json_file, 'r') as json_handle:
            input_vals = json.load(json_handle)

            # Here are the initial heat exchanger values
            self.tube_Length = (input_vals["initial_specs"]["Heat_exchanger"]["L_tube[m]"])
            self.Di = (input_vals["initial_specs"]["Heat_exchanger"]["Di[cm]"])*10**-2 # value should now be in meters
            self.Do = (input_vals["initial_specs"]["Heat_exchanger"]["Do[cm]"])*10**-2 # value should now be in meters
            self.wall_thickness = (input_vals["initial_specs"]["Heat_exchanger"]["wall_thickness[mm]"])*10**-3 # value should now be in meters
            self.wall_material = input_vals["initial_specs"]["Heat_exchanger"]["wall_material"]
            self.nusselt_in = input_vals["initial_specs"]["Heat_exchanger"]["Nusselt_in"]
            self.nusselt_out = input_vals["initial_specs"]["Heat_exchanger"]["Nusselt_out"]
            self.T_in_H2O = fahrenheit_to_kelvin((input_vals["initial_specs"]["Heat_exchanger"]["T_in_H2O[deg_F]"])) # now this value is in kelvin
            self.T_in_wort = fahrenheit_to_kelvin((input_vals["initial_specs"]["Heat_exchanger"]["T_in_wort[deg_F]"])) # now this value is in kelvin
            self.T_out_wort = fahrenheit_to_kelvin((input_vals["initial_specs"]["Heat_exchanger"]["T_out_wort[deg_F]"]))# now this value is in kelvin
            
            # Here are the initial pump values
            self.wort_sum_minor_losses = input_vals["initial_specs"]["pump"]["wort_sum_minor_losses"]
            self.P_in_H2O = psi_to_Pa(input_vals["initial_specs"]["pump"]["P_in_H2O[psi]"]) # this value is now in Pa 
            self.wort_delta_z = feet_to_m(input_vals["initial_specs"]["pump"]["wort_delta_z[ft]"]) # value is now in m. 
            self.density_H2O = input_vals["initial_specs"]["pump"]["density_H2O[kg/m^3]"]
            self.starting_fric_guess = input_vals["initial_specs"]["pump"]["starting_fric_guess"]
            self.H2O_dynamic_viscosity = (input_vals["initial_specs"]["pump"]["H2O_dynamic_viscosity[N-s/m^2]*10^6"])*10**(-6) 

            # Here is an option to turn on the debug (runs debugger function if this input is checked as true)
            self.debug = input_vals["initial_specs"]["debug"]
            self.iteration_stopper = input_vals["initial_specs"]["iteration_stopper"]


    def calc_hydraulic_diameter_annulus(self, Do, thick, Di):
        """This function calculates the hydraulic diameter of the annulus given the following inputs:
        
        Parameters
        ----------
        - Do: Outer tube outer diameter [m] (float)
        - Di: Inner tube outer diameter [m] (float)
        - thick: wall thickness [m] (float)"""
        
        hydraulic_diameter = Do - thick - Di #### this is dependent on if the wall thickness is important here
        return hydraulic_diameter
    
    
    def calc_hydraulic_area_annulus(self, outer_diam_outer_pipe, thick, outer_diam_inner_pipe):
        """This function calculates the hydraulic area of the annulus using the following inputs:
        
        Parameters
        ----------
        - outer_diam_outer_pipe: the outer diameter of the outer pipe that holds the water [m] (float)
        - outer_diam_inner_pipe: the outer diameter of the inner pipe that holds the wort [m] (float)"""

        inner_diam_outer_pipe = outer_diam_outer_pipe - thick
        inner_area_outer_pipe = calc_cyl_cross_area(inner_diam_outer_pipe)
        outer_area_inner_pipe = calc_cyl_cross_area(outer_diam_inner_pipe)
        
        hydraulic_area_annulus = inner_area_outer_pipe - outer_area_inner_pipe
        return hydraulic_area_annulus
    

    def calc_velocity(self, dp, dx, hyd_diam, f, rho):
        """This function returns a tube velocity with the following inputs:
            
            Parameters
            ----------
            - dp: delta pressure [Pa] (float) 
            - dx: tube length [m] (float) 
            - Do: Outer tube flow diameter [m] (float) 
            - f: friction factor (float) 
            - rho: density [kg/m^3] (float)"""
        
        dp_over_dx = dp/dx
        velocity_water = np.sqrt((2*dp_over_dx*hyd_diam)/(f*rho))

        return velocity_water


    def calc_Reynolds(self, rho, v, hyd_diam, mu):
        """This function returns a Reynolds number with the following inputs:
            
            Parameters
            ----------
            - rho: density [kg/m^3] (float)
            - v, velocity [m/s] (float)
            - L, pipe length [m] (float)
            - mu, Dynamic viscosity [N-s/m^2] (float)"""
        
        Reynolds = (rho*v*hyd_diam)/mu
        return Reynolds

    
    def calc_friction_factor(self, Re):
        """This function returns the darcy friction factor given the following inputs:
            
            Parameters
            ----------
            - Re: Reynolds number (float)
            """
        f = (0.79*np.log(Re) - 1.64)**(-2)
        return f 


    # first, we need to calculate the k factor of the 


    def run(self, Do, Di, dp, dx, f, rho):
        """This is a run function that uses all of the functions above.
        
        Parameters
            ----------
            - Do: Outer tube flow diameter [m] (float) 
            - Di: Inner tube outer diameter [m] (float)
            - dp: delta pressure [Pa] (float) 
            - dx: tube length [m] (float) 
            - f: friction factor (float)
            - rho: density [kg/m^3] (float)"""
        
        thickness = self.wall_thickness
        H2O_dyn_visc = self.H2O_dynamic_viscosity
        hyd_diam = self.calc_hydraulic_diameter_annulus(Do, thickness, Di)
        hyd_area = self.calc_hydraulic_area_annulus(Do, thickness, Di)

        water_vel = self.calc_velocity(dp, dx, hyd_diam, f, rho)
        water_Re = self.calc_Reynolds(rho, water_vel, hyd_diam, H2O_dyn_visc)
        water_fric_fact = self.calc_friction_factor(water_Re)
        epsilon = 10 # starter 
        while epsilon >= self.iteration_stopper:
            #  print("iteration!!!!")
             water_vel = self.calc_velocity(dp, dx, hyd_diam, water_fric_fact, rho)
             water_Re = self.calc_Reynolds(rho, water_vel, hyd_diam, H2O_dyn_visc)
             water_fric_fact_new = self.calc_friction_factor(water_Re)
             epsilon = abs((water_fric_fact-water_fric_fact_new)/(water_fric_fact_new))
             water_fric_fact = water_fric_fact_new

        mdot = calc_mdot(rho, water_vel, hyd_area)
        vdot = calc_vdot(rho, mdot)



        if self.debug == True:
            self.debugger(hyd_diam, hyd_area, water_vel, water_Re, water_fric_fact, mdot, vdot)


    def debugger(self, hyd_diam, hyd_area, H2O_vel, H2O_Re, H2O_fric, mdot, vdot):
        """This function prints out every calculation if debug is marked as True in the input file"""

        
        #These are the given values
        print("\n tube_Length [m]:\n",  self.tube_Length, "[m]\n")
        print(" Di [m]:\n",  self.Di, "[m]\n")
        print(" Do [m]:\n",  self.Do, "[m]\n")
        print(" Initial wall thickness [m]:\n",  self.wall_thickness, "[m]\n")
        print(" Initial wall material [string]:\n",  self.wall_material, "\n")
        print(" Nusselt number going in:\n",  self.nusselt_in, "\n")
        print(" Nusselt number going out:\n",  self.nusselt_out, "\n")
        print(" Water temperature in (Tci) [K]:\n",  self.T_in_H2O, "[K]\n")
        print(" Wort temperature in (THi) [K]:\n",  self.T_in_wort, "[K]\n")
        print(" Wort temperature out (THo) [K]:\n",  self.T_out_wort, "[K]\n")
        print(" Sum of minor losses:\n",  self.wort_sum_minor_losses, "\n")
        print(" Water pressure in [Pa]:\n",  self.P_in_H2O, "[Pa]\n")
        print(" Water density [kg/m^3]:\n",  self.density_H2O, "[kg/m^3]\n")
        print(" Wort Delta Z [m]:\n",  self.wort_delta_z, "[m]\n")
        print(" Starting darcy friction factor, f:\n",  self.starting_fric_guess, "\n")
        print(" Water Dynamic Viscosity, f:\n",  self.H2O_dynamic_viscosity, "[m]\n")

        #Here is where we start making calculations
        print(" Hydraulic Diameter of area between tubes (space carrying water):\n",  hyd_diam, "[m]\n")
        print(" Hydraulic cross sectional area between tubes (space carrying water):\n",  hyd_area, "[m^2]\n")
        print(" Water velocity [m/s]:\n",  H2O_vel, "[m]\n")
        print(" Water Reynolds :\n",  H2O_Re, "\n")
        print(" Water Darcy Friction Factor :\n",  H2O_fric, "\n")
        print(" Mass flow rate of the water going through the pipe [kg/s]:\n",  mdot, "[kg/s]\n")
        print(" Volume flow rate of the water going through the pipe [m^3/s]:\n",  vdot, "[m^3/s]\n")
        

        
if __name__ == "__main__":
    choice = pump("pump_selection.json")
    choice.run(choice.Do, choice.Di, choice.P_in_H2O, choice.tube_Length, choice.starting_fric_guess, choice.density_H2O) 
    
