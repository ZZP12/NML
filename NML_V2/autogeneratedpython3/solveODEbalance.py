import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import math

from dataDictionary import generate_model_parameters_dictionary
from kinetics import calculate_kinetics
from ODEbalance import simulation_ODEs

time_start = 0.0
time_stop = 10.0
time_num = int(100 * (time_stop - time_start))
t = np.linspace(time_start, time_stop, time_num)

# initial species concentration
y0 = [
	0.001,  # 0  BIOMASS
	0.001,  # 1  met_A_e
	0.001,  # 2  met_B_e
	0.001,  # 3  met_A_c
	0.001,  # 4  met_B_c
	0.001,  # 5  p_C_c
	0.001,  # 6  p_E_c
	0.0     # 7  mRNA_E_c
]

dataDictionary = generate_model_parameters_dictionary()
all_species_reversed_dict = dataDictionary["all_species_reversed_dict"]

y = integrate.odeint(simulation_ODEs, y0, t, args=(dataDictionary,))

# plot results
col = len(y0)
subplt_col = math.ceil(col/4)
for i in range(col):
	plt.subplot(4, subplt_col, i+1)
	plt.plot(t, y[:, i])
	plt.title(all_species_reversed_dict[i])
plt.show()