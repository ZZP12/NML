include("./include.jl")

time_start = 0.0
time_step_size = 0.1
time_stop = 10.0
t = collect(time_start:time_step_size:time_stop)

# initial species concentration
y0 = [
	0.001,  # BIOMASS
	0.001,  # met_A_e
	0.001,  # met_B_e
	0.001,  # met_A_c
	0.001,  # met_B_c
	0.001,  # p_C_c
	0.001,  # p_E_c
	0.0     # mRNA_E_c
]

dataDictionary = generate_model_parameters_dictionary()
all_species_reversed_dict = dataDictionary["all_species_reversed_dict"]

f(t, y) = simulation_ODEs(t, y, dataDictionary)
t, y = ode23s(f, y0, t)

# data transfer
row = length(t)
col = length(y0)
Y = zeros(row, col)
foreach(x->(Y[x, :] = y[x]), collect(1:row))

# plot results
subplt_col = round(Int, col/4)
for i = 1:col
	subplot(4, subplt_col, i)
	plot(t, Y[:,i])
	title(all_species_reversed_dict[i])
end