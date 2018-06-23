function generate_model_parameters_dictionary()

	# load stoichiometry
	stoichiometric_matrix = readdlm("./stoichiometry.dat")

	# all species dictionary
	all_species_dict = Dict{String, Int64}()  # String --> Int
	all_species_dict["met_A_e"] = 1
	all_species_dict["met_B_e"] = 2
	all_species_dict["met_A_c"] = 3
	all_species_dict["met_B_c"] = 4
	all_species_dict["p_C_c"] = 5
	all_species_dict["p_E_c"] = 6
	all_species_dict["mRNA_E_c"] = 7

	# all species dictionary_reversed
	all_species_reversed_dict = Dict{Int64, String}()  # Int --> String
	for (key,val) in all_species_dict
		all_species_reversed_dict[val] = key
	end

	# all rnx species array
	rnx_species_array::Array{String,1} = [
		"met_A_e",
		"met_B_e",
		"met_A_c",
		"met_B_c",
		"p_C_c",
		"p_E_c"
	]

	# all mRNA species array
	mRNA_species_array::Array{String,1} = [
		"mRNA_E_c"
	]

	# all protein species array
	protein_species_array::Array{String,1} = [
		"p_E_c"
	]

	# signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)
	kcat_signaling = 1e-3*ones(5)  # kcat[#reaction]: reaction name
	kcat_signaling[1] = 1.0  # kcat: 1.0*met_A_e<uptake:>1.0*met_A_c
	kcat_signaling[2] = 1.0  # kcat: 1.0*met_A_c<catalyze:1.0*p_E_c>1.0*met_B_c
	kcat_signaling[3] = 1.0  # kcat: 1.0*met_B_c<secrete:>1.0*met_B_e
	kcat_signaling[4] = 1.0  # kcat: 1.0*p_C_c<react:>
	kcat_signaling[5] = 1.0  # kcat: <react:>1.0*p_C_c

	# Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate
	MonodAffinityConstantDict = Dict{String, Float64}()
	MonodAffinityConstantDict["MonodK_1_met_A_e"] = 1.0
	MonodAffinityConstantDict["MonodK_2_met_A_c"] = 1.0
	MonodAffinityConstantDict["MonodK_3_met_B_c"] = 1.0
	MonodAffinityConstantDict["MonodK_4_p_C_c"] = 1.0


	# W_value in transcription control term, nomenclature: W_targetmRNA_actor
	W_value_dict = Dict{String, Float64}()
	W_value_dict["W_mRNA_E_c_p_C_c"] = 1.0

	# disassociation constants in transfer function, nomenclature: KD_speciesName
	transferFunctionDisassociationConstantDict = Dict{String, Float64}()
	transferFunctionDisassociationConstantDict["KD_p_C_c"] = 1.0

	# transcription specific correction factor
	transcriptionSpecificCorrectionFactorDict = Dict{String, Float64}()
	transcriptionSpecificCorrectionFactorDict["mRNA_E_c"] = 1.0  # mRNA_E_c

	# background gene expression control term
	backgroundGeneExpressionControlTermDict = Dict{String, Float64}()
	backgroundGeneExpressionControlTermDict["mRNA_E_c"] = 0.1  # mRNA_E_c

	# translation specific correction factor
	translationSpecificCorrectionFactor = Dict{String, Float64}()
	translationSpecificCorrectionFactor["p_E_c"] = 1.0  # p_E_c

	# cooperativity number in transfer function
	cooperativity = 1


	# ------------------------------------------------------------------------------------------#
	# constants (from bionumbers)       units
	# ------------------------------------------------------------------------------------------#
	cell_diameter = 1.1                 # mum
	number_of_rnapII = 4600            	# copies/cells
	number_of_ribosome = 50000         	# copies/cells
	mRNA_half_life_TF = 0.083           # hrs
	protein_half_life = 70              # hrs
	doubling_time_cell = 0.33           # hrs
	max_translation_rate = 16.5         # aa/sec
	max_transcription_rate = 60.0       # nt/sec
	average_transcript_length = 1200   	# nt
	average_protein_length = 400       	# aa
	fraction_nucleus = 0.0             	# dimensionless
	av_number = 6.02e23                 # number/mol
	avg_gene_number = 2                 # number of copies of a gene
	polysome_number = 4									# number of ribsomoses per transcript
	# ------------------------------------------------------------------------------------------#
	#
	# ------------------------------------------------------------------------------------------#
	# Calculate constants using bionumber values
	# ------------------------------------------------------------------------------------------#
	# Calculate the volume (convert to L)
	V = ((1-fraction_nucleus)*(1/6)*(3.14159)*(cell_diameter)^3)*(1e-15)
	
	# Calculate the rnapII_concentration and ribosome_concentration
	rnapII_concentration = number_of_rnapII*(1/av_number)*(1/V)*1e9                   # nM
	ribosome_concentration = number_of_ribosome*(1/av_number)*(1/V)*1e9               # nM
	
	# degrdation rate constants -
	degradation_constant_mRNA = -(1/mRNA_half_life_TF)*log(0.5)                       # hr^-1
	degradation_constant_protein = -(1/protein_half_life)*log(0.5)                    # hr^-1
	
	# kcats for transcription and translation -
	kcat_transcription = max_transcription_rate*(3600/average_transcript_length)      # hr^-1
	kcat_translation = polysome_number*max_translation_rate*(3600/average_protein_length)             # hr^-1
	
	# Maximum specific growth rate -
	maximum_specific_growth_rate = (1/doubling_time_cell)*log(2)                      # hr^-1
	
	# What is the average gene concentration -
	avg_gene_concentration = avg_gene_number*(1/av_number)*(1/V)*1e9                  # nM
	
	# How fast do my cells die?
	death_rate_constant = 0.05*maximum_specific_growth_rate                            # hr^-1
	
	# Saturation constants for translation and trascription -
	saturation_transcription = 4600*(1/av_number)*(1/V)*1e9                           # nM
	saturation_translation = 150000*(1/av_number)*(1/V)*1e9                           # nM
	# -------------------------------------------------------------------------------------------#


	####################################
	# put all stuff in a Dictionary
	dataDictionary = Dict{String, Any}()
	dataDictionary["stoichiometric_matrix"] = stoichiometric_matrix
	dataDictionary["all_species_dict"] = all_species_dict
	dataDictionary["all_species_reversed_dict"] = all_species_reversed_dict
	dataDictionary["rnx_species_array"] = rnx_species_array
	dataDictionary["mRNA_species_array"] = mRNA_species_array
	dataDictionary["protein_species_array"] = protein_species_array
	dataDictionary["kcat_signaling"] = kcat_signaling
	dataDictionary["MonodAffinityConstantDict"] = MonodAffinityConstantDict
	dataDictionary["W_value_dict"] = W_value_dict
	dataDictionary["transferFunctionDisassociationConstantDict"] = transferFunctionDisassociationConstantDict
	dataDictionary["transcriptionSpecificCorrectionFactorDict"] = transcriptionSpecificCorrectionFactorDict
	dataDictionary["backgroundGeneExpressionControlTermDict"] = backgroundGeneExpressionControlTermDict
	dataDictionary["translationSpecificCorrectionFactor"] = translationSpecificCorrectionFactor
	dataDictionary["cooperativity"] = cooperativity
	dataDictionary["RNAPConcentration"] = rnapII_concentration
	dataDictionary["RIBOConcentration"] = ribosome_concentration
	dataDictionary["degradationConstantmRNA"] = degradation_constant_mRNA
	dataDictionary["degradationConstantProtein"] = degradation_constant_protein
	dataDictionary["kcatTranscription"] = kcat_transcription
	dataDictionary["kcatTranslation"] = kcat_translation
	dataDictionary["avgGeneConcentration"] = avg_gene_concentration
	dataDictionary["transcriptionSaturationConstant"] = saturation_transcription
	dataDictionary["translationSaturationConstant"] = saturation_translation
	dataDictionary["specificGrowthRate"] = maximum_specific_growth_rate - death_rate_constant

	return dataDictionary

end