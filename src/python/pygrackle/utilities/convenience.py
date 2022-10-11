########################################################################
#
# Python wrapper convenience functions
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np
import sys

from pygrackle.fluid_container import \
    FluidContainer

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, mass_electron_cgs,boltzmann_constant_cgs

def check_convergence(fc1, fc2, fields=None, tol=0.01):
    "Check for fields to be different by less than tol."

    if fields is None:
        fields = ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                  "H2I", "H2II", "DI", "DII", "HDI", "de"]
    max_field = None
    max_val = 0.0
    for field in fields:
        if field not in fc1 or field not in fc2:
            continue
        convergence = np.max(np.abs(fc1[field] - fc2[field]) / fc1[field])
        if convergence > max_val:
            max_val = convergence
            max_field = field
    if np.any(max_val > tol):
        sys.stderr.write("max change - %5s: %.10e." % (max_field, max_val))
        return False
    return True

def setup_fluid_container(my_chemistry,
                          density=mass_hydrogen_cgs,
                          temperature=None,
                          hydrogen_mass_fraction=0.76,
                          metal_mass_fraction=0.02041,
                          d_to_h_ratio=3.4e-5,
                          converge=False, tolerance=0.01,
                          max_iterations=1):
    """
    Initialize a fluid container with a constant density and smoothly
    increasing temperature from 10 K to 1e9 K.  Optionally, iterate the
    chemistry solver until the species fractions converge.  Return
    the fluid container.

    The input are expected to be in CGS.
    """

    rval = my_chemistry.initialize()
    if rval == 0:
        raise RuntimeError("Failed to initialize chemistry_data.")

    tiny_number = 1e-13
    if temperature is None:
        n_points = 200
        temperature = np.logspace(4, 9, n_points)
    else:
        n_points = temperature.size
    fc = FluidContainer(my_chemistry, n_points)
    fc["density"][:] = density / my_chemistry.density_units
    if my_chemistry.primordial_chemistry > 0:
        fc["HI"][:] = hydrogen_mass_fraction * fc["density"]
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - hydrogen_mass_fraction) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        fc["de"][:] = tiny_number * fc["density"]#fc["de"][:] - fc["HM"] + fc["H2II"] / 2.0    
    if my_chemistry.primordial_chemistry > 1:
        """    
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
        fc["H2I"][:] = hydrogen_mass_fraction * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        """
        fc["HI"][:] = 0.5*hydrogen_mass_fraction * fc["density"] #tiny_number * fc["density"]
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
        fc["H2I"][:] = 0.5*hydrogen_mass_fraction * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - hydrogen_mass_fraction) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        fc["de"][:] = tiny_number * fc["density"]#fc["de"][:] - fc["HM"] + fc["H2II"] / 2.0    
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * d_to_h_ratio * \
        fc["density"]
        #(fc["HI"][:]+fc["HII"][:]+fc["HM"][:]+fc["H2I"][:]+fc["H2II"][:])
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    fc["metal"][:] = metal_mass_fraction * fc["density"]
    fc["dust"][:] = 0.009387 * fc["density"]


    if my_chemistry.primordial_chemistry > 0:
        fc["de"][:] =  fc["de"][:]*(mass_hydrogen_cgs/mass_electron_cgs)

    print("mu 1 = ",fc["mu"])
    print("gamma 1 = ",fc["gamma"])
    print("temperature 1 = ",fc["temperature"])
    fc["temperature"] = temperature
    fc.calculate_gamma()
    fc.calculate_mean_molecular_weight()
    fc["energy"] = temperature / \
    fc.chemistry_data.temperature_units / \
    fc["mu"] / (fc["gamma"] - 1.0)
    
    print("mu B = ",fc["mu"])
    print("gamma B = ",fc["gamma"])
    print("temperature B = ",temperature)
    
    fc.calculate_gamma()

    print("mu C = ",fc["mu"])
    print("gamma C = ",fc["gamma"])
    print("temperature C = ",temperature)

    fc["energy"] = temperature / \
    fc.chemistry_data.temperature_units / \
    fc["mu"] / (fc["gamma"] - 1.0)

    #if my_chemistry.primordial_chemistry > 0:
    fc.calculate_mean_molecular_weight()
    fc.calculate_gamma()
    #print("mu 2b = ",fc["mu"])
    #print("gamma 2b = ",fc["gamma"])        
    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    if my_chemistry.primordial_chemistry > 0:
        fc["specific_heating_rate"][:] = 0.
        fc["volumetric_heating_rate"][:] = 0.    

    fc_last = fc.copy()

    my_time = 0.0
    i = 0

    if my_chemistry.primordial_chemistry > 0:
        fc["de"][:] =  fc["de"][:]*(mass_electron_cgs/mass_hydrogen_cgs)    
    while converge and i < max_iterations:
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_hydrogen_cgs/mass_electron_cgs)        
        fc.calculate_cooling_time()
        #if my_chemistry.primordial_chemistry > 0:
        #    fc.calculate_gamma()
        #fc.calculate_pressure()
        dt = 0.1 * np.abs(fc["cooling_time"]).min()
        #dt = 2.230e12/my_chemistry.time_units
        #sys.stderr.write("t: %.3f Myr, dt: %.3e Myr, " % \
        #                 ((my_time * my_chemistry.time_units / sec_per_Myr),
        #                  (dt * my_chemistry.time_units / sec_per_Myr)))
        print("t: %.3e s, dt: %.3e s, " % \
                         ((my_time * my_chemistry.time_units),
                          (dt * my_chemistry.time_units)))        
        for field in ["HI", "HII", "HM", "HeI", "HeII", "HeIII",
                      "H2I", "H2II", "DI", "DII", "HDI", "de"]:
            if field in fc:
                fc_last[field] = np.copy(fc[field])
        print("\n", "Before :", "\n")
        print("velocity_units :",my_chemistry.velocity_units)
        print("velocity_units**2 :",my_chemistry.velocity_units**2)
        print("Temperature : ", fc["temperature"][-1])
        print("energy : ", fc["energy"][-1]*my_chemistry.energy_units)
        print("density : ", fc["density"][-1]*my_chemistry.density_units)
        
        print("pressure : ", fc["pressure"][-1]*my_chemistry.pressure_units)
        print("mu : ", fc["mu"][-1])
        print("gamma : ", fc["gamma"][-1])
        print("cooling_time : ", fc["cooling_time"][-1]*my_chemistry.time_units)
        print("dt : ", dt*my_chemistry.time_units)
        if my_chemistry.primordial_chemistry > 0:
            for field in ["HI", "HII",  "HeI", "HeII", "HeIII",
                         "de", "metal", "dust"]:
                print(field," : ",fc[field][-1]* my_chemistry.density_units)
        if my_chemistry.primordial_chemistry > 1:                           
            for field in ["HM", "H2I", "H2II"]:
                print(field," : ",fc[field][-1]* my_chemistry.density_units)             
        if my_chemistry.primordial_chemistry > 2:
            for field in ["DI","DII","HDI"]:                        
                print(field," : ",fc[field][-1]* my_chemistry.density_units)
        
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_electron_cgs/mass_hydrogen_cgs)             
        density_proper = fc["density"] / \
        (my_chemistry.a_units *
         my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)
        cooling_rate = fc.chemistry_data.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper
        print("cooling_rate : ", cooling_rate[-1])
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_hydrogen_cgs/mass_electron_cgs)          
        print("energy 1 : ", fc["energy"][-1]*my_chemistry.energy_units)
        print("cooling time 1 : ", fc["cooling_time"][-1]*my_chemistry.time_units)
        print("mu 1 : ", fc["mu"][-1])
        fc.solve_chemistry(dt)
        fc.calculate_cooling_time()
        fc.calculate_gamma()
        fc.calculate_pressure()
        fc.calculate_temperature()
        if my_chemistry.primordial_chemistry == 0:
            fc["temperature"][:]=fc["pressure"][:]*my_chemistry.pressure_units*\
                mass_hydrogen_cgs/(boltzmann_constant_cgs*\
                   fc["density"][:]*my_chemistry.density_units)
        print("\n", "After :", "\n")
        print("Temperature : ", fc["temperature"][-1])
        print("energy : ", fc["energy"][-1]*my_chemistry.energy_units)
        print("density : ", fc["density"][-1]*my_chemistry.density_units)
        print("pressure : ", fc["pressure"][-1]*my_chemistry.pressure_units)
        print("mu : ", fc["mu"][-1])
        print("gamma : ", fc["gamma"][-1])
        print("cooling_time : ", fc["cooling_time"][-1]*my_chemistry.time_units)
        print("dt : ", dt*my_chemistry.time_units)
        if my_chemistry.primordial_chemistry > 0:
            for field in ["HI", "HII",  "HeI", "HeII", "HeIII",
                         "de", "metal", "dust"]:
                print(field," : ",fc[field][-1]* my_chemistry.density_units)
        if my_chemistry.primordial_chemistry > 1:                           
            for field in ["HM", "H2I", "H2II"]:
                print(field," : ",fc[field][-1]* my_chemistry.density_units)             
        if my_chemistry.primordial_chemistry > 2:
            for field in ["DI","DII","HDI"]:                        
                print(field," : ",fc[field][-1]* my_chemistry.density_units)
                        
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_electron_cgs/mass_hydrogen_cgs)                
        density_proper = fc["density"] / \
        (my_chemistry.a_units *
         my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)
        cooling_rate_print = fc.chemistry_data.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper
        print("cooling_rate 1 : ", cooling_rate_print[-1])
        print("energy 2 : ", fc["energy"][-1]*my_chemistry.energy_units)
        print("cooling time 2 : ", fc["cooling_time"][-1]*my_chemistry.time_units)
        print("mu 2 : ", fc["mu"][-1])
        #print("mu 3 = ",fc["mu"])
        #print("gamma 3 = ",fc["gamma"])
        #print("temperature 3 = ",fc["temperature"])
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_hydrogen_cgs/mass_electron_cgs)         
        fc.calculate_mean_molecular_weight()
        if my_chemistry.primordial_chemistry > 0:
            fc["de"][:] =  fc["de"][:]*(mass_electron_cgs/mass_hydrogen_cgs)            
        #print("mu 4 = ",fc["mu"])
        #print("gamma 4 = ",fc["gamma"])
        #print("temperature 4 = ",fc["temperature"])
        fc["energy"] = temperature / \
            fc.chemistry_data.temperature_units / fc["mu"] / \
            (my_chemistry.Gamma - 1.0)
        my_time += dt
        print("t: %.3e s, dt: %.3e s, " % \
                         ((my_time * my_chemistry.time_units),
                          (dt * my_chemistry.time_units)))            
        converged = check_convergence(fc, fc_last, tol=tolerance)
        if converged:
            sys.stderr.write("\n")
            break
        sys.stderr.write("\r")

        i += 1

    """if i >= max_iterations:
        sys.stderr.write("ERROR: solver did not converge in %d iterations.\n" %
                         max_iterations)
        return None
    """

    print("final time = %.8e s" % (my_time*my_chemistry.time_units))
    return fc
