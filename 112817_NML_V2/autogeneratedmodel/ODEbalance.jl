# ----------------------------------------------------------------------------------- #
# Copyright (c) 2017 Varnerlab
# Robert Frederick Smith School of Chemical and Biomolecular Engineering
# Cornell University, Ithaca NY 14850
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ----------------------------------------------------------------------------------- #


# set up ODE & DAE
function simulation_ODEs(t,y,dataDictionary)
	# get the derivatives

	# correct for negatives -
	idx_small = find(y.<0)
	y[idx_small] = 0.0

	# load data from data dictionary
	rnx_species::Array{String,1} = dataDictionary["rnx_species_array"]
	mRNA_species::Array{String,1} = dataDictionary["mRNA_species_array"]
	protein_species::Array{String,1} = dataDictionary["protein_species_array"]
	stoichiometry = dataDictionary["stoichiometric_matrix"]
	all_species_dict::Dict{String, Int} = dataDictionary["all_species_dict"]  # String --> number
	specific_growth_rate = dataDictionary["specificGrowthRate"]
	transcription_correction = dataDictionary["transcriptionSpecificCorrectionFactorDict"]
	translation_correction = dataDictionary["translationSpecificCorrectionFactor"]
	degradation_constant_mRNA = dataDictionary["degradationConstantmRNA"]
	degradation_constant_protein = dataDictionary["degradationConstantProtein"]

	# kinetics
	(rnx_rate, transcription_rate, translation_rate) = calculate_kinetics(y, dataDictionary)

	# calculate the derivatives
	dydt = zeros(length(all_species_dict))

	# signaling reaction network
	X2 = zeros(length(rnx_species)-2)  # species inside the cell
	for (id, st) in enumerate(rnx_species[2+1:end])  # get X2 value from y
		X2[id] = y[all_species_dict[st]]
	end
	dX1 = stoichiometry[1:2, :]*rnx_rate  # extracellular metabolites
	dX2 = stoichiometry[2+1:end, :]*rnx_rate - specific_growth_rate*X2
	dX = [dX1; dX2]
	for (id, st) in enumerate(rnx_species)  # return derivatives back to dydt
		dydt[all_species_dict[st]] = dX[id]
	end

	# TXTL network
	for i = 1:length(mRNA_species)  # seems easy to add additional terms to account for the signaling effect
		mRNA_id = all_species_dict[mRNA_species[i]]
		p_id = all_species_dict[protein_species[i]]
		dydt[mRNA_id] = transcription_rate[mRNA_species[i]] - (specific_growth_rate + transcription_correction[mRNA_species[i]]*degradation_constant_mRNA)*y[mRNA_id]
		dydt[p_id] = translation_rate[protein_species[i]] - (specific_growth_rate+ translation_correction[protein_species[i]]*degradation_constant_protein)*y[p_id]
	end

	return dydt  # derivative
end


function simulation_DAEs(t, y, yp, r, dataDictionary)
	dydt = simulation_ODEs(t, y, dataDictionary)
	for i = 1:length(dydt)
		r[i] = dydt[i] - yp[i]
	end
end