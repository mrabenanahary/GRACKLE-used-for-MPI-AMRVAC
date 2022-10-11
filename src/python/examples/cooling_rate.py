########################################################################
#
# Cooling rate example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import numpy as np
import os
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

if __name__ == "__main__":
    current_redshift = 0.

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    if 'PRIMORDIAL_CHEM' in os.environ:
        my_chemistry.primordial_chemistry = int(os.environ['PRIMORDIAL_CHEM'])
    else:
        my_chemistry.primordial_chemistry = 1
    number_of_iter = 1
    my_chemistry.dust_chemistry = 1
    my_chemistry.h2_on_dust = 0
    my_chemistry.metal_cooling = 1
    my_chemistry.use_dust_density_field = 1
    my_chemistry.UVbackground = 0
    my_chemistry.self_shielding_method = 0
    my_chemistry.H2_self_shielding = 0
    my_dir = os.path.dirname(os.path.abspath(__file__))
    grackle_data_file = bytearray(os.path.join(
        my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012_high_density.h5"), 'utf-8')
    my_chemistry.grackle_data_file = grackle_data_file

    my_chemistry.use_specific_heating_rate = 1
    my_chemistry.use_volumetric_heating_rate = 1

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1.0 / (1.0 + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Gyr in s
    print("density_units", my_chemistry.density_units)
    print("length_units", my_chemistry.length_units)
    print("time_units", my_chemistry.time_units)
    my_chemistry.set_velocity_units()

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    #temperature = np.array([1.00000000e+09])
    temperature = np.logspace(1, 9, 200)

    temperature_zero = np.copy(temperature)

    make_converge = True
    
    fc = setup_fluid_container(my_chemistry,
                               density=mass_hydrogen_cgs*1.0e0,
                               temperature=temperature,
                               converge=make_converge,
                               max_iterations=number_of_iter)

    if my_chemistry.primordial_chemistry > 0:
        fc["specific_heating_rate"][:] = 0.
        fc["volumetric_heating_rate"][:] = 0.

    
    if(make_converge==False):
        fc.calculate_temperature()
        fc.calculate_pressure()
        fc.calculate_cooling_time()

    density_proper = fc["density"] / \
        (my_chemistry.a_units *
         my_chemistry.a_value)**(3*my_chemistry.comoving_coordinates)
    cooling_rate = fc.chemistry_data.cooling_units * fc["energy"] / \
        np.abs(fc["cooling_time"]) / density_proper
    print("cooling_rate 2: ", cooling_rate[-1])
    print("energy 3 : ", fc["energy"][-1]*my_chemistry.energy_units)
    print("cooling time 3 : ", fc["cooling_time"][-1]*my_chemistry.time_units)
    print("mu 3 : ", fc["mu"][-1])

    #print("mu =")
    #print(fc["mu"])

    print("temperautre zero = ",temperature_zero)
    data = {}
    t_sort = np.argsort(fc["temperature"])
    for field in fc.density_fields:
        data[field] = yt.YTArray(fc[field][t_sort] *
                                 my_chemistry.density_units, "g/cm**3")
    data["energy"] = yt.YTArray(
        fc["energy"][t_sort] * my_chemistry.energy_units, "erg/g")
    data["temperature"] = yt.YTArray(fc["temperature"][t_sort], "K")
    data["temperature_zero"] = yt.YTArray(temperature_zero[t_sort], "K")
    data["pressure"] = yt.YTArray(
        fc["pressure"][t_sort] * my_chemistry.pressure_units, "dyne/cm**2")
    data["cooling_time"] = yt.YTArray(fc["cooling_time"][t_sort]*my_chemistry.time_units, "s")
    data["cooling_rate"] = yt.YTArray(cooling_rate[t_sort], "erg*cm**3/s")
    data["cooling_rate_zero"] = yt.YTArray(cooling_rate, "erg*cm**3/s")
    data["density_proper"] = yt.YTArray(density_proper[t_sort], "g/cm**3")
    data["mu"] = yt.YTArray(fc["mu"][t_sort],"")
    data["gamma"] = yt.YTArray(fc["gamma"][t_sort],"")

    pyplot.loglog(data["temperature"], data["cooling_rate"],
                  color="black")
    pyplot.xlabel('T [K]')
    pyplot.ylabel('$\\Lambda$ [erg s$^{-1}$ cm$^{3}$]')

    # save data arrays as a yt dataset
    if 'PRIMORDIAL_CHEM' in os.environ:
        ds_name = 'cooling_rate.pc%s.h5' % os.environ['PRIMORDIAL_CHEM']
        im_name = 'cooling_rate.pc%s.png' % os.environ['PRIMORDIAL_CHEM']
    else:
        ds_name = 'cooling_rate.h5'
        im_name = 'cooling_rate.png'
    pyplot.tight_layout()
    pyplot.savefig(im_name)
    yt.save_as_dataset({}, ds_name, data)
