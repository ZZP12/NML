function [rnx_rate_vector, TX_rate_vector, TL_rate_vector] = calculate_kinetics(X, data_dictionary)

	% generate kinetics for signaling, assume saturation kinetics
	% load Monod affinity constants
	MonodK = data_dictionary('MonodAffinityConstantDict');  % String-->Float
	% load reaction kinetic constants
	kcat = data_dictionary('kcat_signaling');  % in order of rnx

	% write reaction rate equations
	rnx_rate_vector = zeros(5, 1);
	rnx_rate_vector(1) = kcat(1)*X(2)/(MonodK('MonodK~rnx1~met_A_e')+X(2));
	rnx_rate_vector(2) = kcat(2)*X(7)*X(4)/(MonodK('MonodK~rnx2~met_A_c')+X(4));
	rnx_rate_vector(3) = kcat(3)*X(5)/(MonodK('MonodK~rnx3~met_B_c')+X(5));
	rnx_rate_vector(4) = kcat(4)*X(6)/(MonodK('MonodK~rnx4~p_C_c')+X(6));
	rnx_rate_vector(5) = kcat(5);

	% generate kinetics & control terms for TXTL;
	% load data from data dictionary;
	W_value_dict = data_dictionary('W_value_dict');  % String-->Float
	coop = data_dictionary('cooperativity');
	disassociation_constant_dict = data_dictionary('transferFunctionDisassociationConstantDict');
	background_control_term = data_dictionary('backgroundGeneExpressionControlTermDict');
	kcatTX = data_dictionary('kcatTranscription');
	RNAPconc = data_dictionary('RNAPConcentration');
	geneConc = data_dictionary('avgGeneConcentration');
	saturationTX = data_dictionary('transcriptionSaturationConstant');
	kcatTL = data_dictionary('kcatTranslation');
	RIBOconc = data_dictionary('RIBOConcentration');
	saturationTL = data_dictionary('translationSaturationConstant');

	% write control terms, transcription rate equations, and translation rate equations
	TX_control_term = containers.Map();  % control term
	TX_rate_vector = containers.Map();  % transcription rate
	TL_rate_vector = containers.Map();  % translation rate
	% mRNA_E_c and p_E_c
	mRNA_E_c_p_C_c = W_value_dict('W~p_C_c~mRNA_E_c') * X(6)^coop/((disassociation_constant_dict('KD~p_C_c'))^coop+X(6)^coop);
	up_action_on_mRNA_E_c = mRNA_E_c_p_C_c ;
	TX_control_term('mRNA_E_c') = (background_control_term('mRNA_E_c') + up_action_on_mRNA_E_c)/(1 + background_control_term('mRNA_E_c') + up_action_on_mRNA_E_c);
	TX_rate_vector('mRNA_E_c') = TX_control_term('mRNA_E_c') * kcatTX * RNAPconc * geneConc / (saturationTX + geneConc);
	TL_rate_vector('p_E_c') = kcatTL * RIBOconc * X(8) / (saturationTL + X(8));

end