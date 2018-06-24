function calculate_kinetics(X, data_dictionary)

	# load all species dictionary
	all_species_dict = data_dictionary["all_species_dict"]  # String-->Int

	# generate kinetics for signaling, assume saturation kinetics
	# load Monod affinity constants
	MonodK = data_dictionary["MonodAffinityConstantDict"]  # String-->Float
	# load reaction kinetic constants
	kcat = data_dictionary["kcat_signaling"]  # in order of rnx

	# write reaction rate equations
	rnx_rate_vector = zeros(7)
	rnx_rate_vector[1] = kcat[1]*X[2]/(MonodK["MonodK~rnx1~m_A_e"]+X[2])
	rnx_rate_vector[2] = kcat[2]*X[5]/(MonodK["MonodK~rnx2~m_A_c"]+X[5])
	rnx_rate_vector[3] = kcat[3]*X[6]/(MonodK["MonodK~rnx3~m_B_c"]+X[6])
	rnx_rate_vector[4] = kcat[4]*X[5]/(MonodK["MonodK~rnx4~m_A_c"]+X[5])
	rnx_rate_vector[5] = kcat[5]*X[7]/(MonodK["MonodK~rnx5~m_C_c"]+X[7])
	rnx_rate_vector[6] = kcat[6]*X[6]/(MonodK["MonodK~rnx6~m_B_c"]+X[6])
	rnx_rate_vector[7] = kcat[7]*X[7]/(MonodK["MonodK~rnx7~m_C_c"]+X[7])

	# generate kinetics & control terms for TXTL
	# load data from data dictionary
	W_value_dict = data_dictionary["W_value_dict"]  #String-->Float
	coop = data_dictionary["cooperativity"]
	disassociation_constant_dict = data_dictionary["transferFunctionDisassociationConstantDict"]
	background_control_term = data_dictionary["backgroundGeneExpressionControlTermDict"]
	kcatTX = data_dictionary["kcatTranscription"]
	RNAPconc = data_dictionary["RNAPConcentration"]
	geneConc = data_dictionary["avgGeneConcentration"]
	saturationTX = data_dictionary["transcriptionSaturationConstant"]
	kcatTL = data_dictionary["kcatTranslation"]
	RIBOconc = data_dictionary["RIBOConcentration"]
	saturationTL = data_dictionary["translationSaturationConstant"]

	# write control terms, transcription rate equations, and translation rate equations
	TX_control_term = Dict{String, Float64}()  # control term
	TX_rate_vector = Dict{String, Float64}()  # transcription rate
	TL_rate_vector = Dict{String, Float64}()  # translation rate

	return (rnx_rate_vector, TX_rate_vector, TL_rate_vector)
end