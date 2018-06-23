clc;
time_start = 0.0;
time_step_size = 0.1;
time_stop = 10.0;
t = time_start:time_step_size:time_stop;

% initial species concentration
y0 = [
	0.001, ...  % BIOMASS
	0.001, ...  % met_A_e
	0.001, ...  % met_B_e
	0.001, ...  % met_A_c
	0.001, ...  % met_B_c
	0.001, ...  % p_C_c
	0.001, ...  % p_E_c
	0.0 ...    % mRNA_E_c
];

dataDictionary = generate_dataDictionary();
all_species_reversed_dict = dataDictionary('all_species_reversed_dict');

[t, y] = ode23s(@(t,y) ODEbalance(t,y,dataDictionary), t, y0);

% plot results
subplt_col = round(length(y0)/4);
for i = 1:length(y0)
	subplot(4, subplt_col, i)
	plot(t, y(:, i))
	title(char(all_species_reversed_dict(i)), 'Interpreter', 'none')
end