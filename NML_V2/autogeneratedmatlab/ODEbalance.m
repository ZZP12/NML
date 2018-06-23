
% set up ODE, get the derivatives
function dydt = ODEbalance(t,y,dataDictionary)
  % correct for negatives
  y(find(y < 0)) = 0.0;

  % load data from data dictionary
  rnx_species = dataDictionary('rnx_species_array');
  mRNA_species = dataDictionary('mRNA_species_array');
  protein_species = dataDictionary('protein_species_array');
  stoichiometry = dataDictionary('stoichiometric_matrix');
  all_species_dict = dataDictionary('all_species_dict');  % String --> number
  specific_growth_rate = dataDictionary('specificGrowthRate');
  transcription_correction = dataDictionary('transcriptionSpecificCorrectionFactorDict');
  translation_correction = dataDictionary('translationSpecificCorrectionFactor');
  degradation_constant_mRNA = dataDictionary('degradationConstantmRNA');
  degradation_constant_protein = dataDictionary('degradationConstantProtein');

  % kinetics
  [rnx_rate, transcription_rate, translation_rate] = calculate_kinetics(y, dataDictionary);

  % calculate the derivatives
  dydt = zeros(1, length(all_species_dict));

  % signaling reaction network
  X2 = zeros(1, (length(rnx_species)-2));  % species inside the cell
  for i = 3:length(rnx_species)  % get X2 value from y
    X2(i - 2) = y(all_species_dict(char(rnx_species(i))));
  end
  dX1 = stoichiometry(1:2, :) * rnx_rate;  % extracellular metabolites
  dX2 = stoichiometry(2+1:end, :)*rnx_rate - transpose(specific_growth_rate*X2);  % INTRA
  dX = [dX1; dX2];
  for i = 1:length(rnx_species)  % return derivatives back to dydt
    dydt(all_species_dict(char(rnx_species(i)))) = dX(i);
  end
  dydt(1) = specific_growth_rate * y(1);  % for Biomass

  % TXTL network
  for i = 1:length(mRNA_species)
  % seems easy to add additional terms to account for the signaling effect
    mRNA_id = all_species_dict(char(mRNA_species(i)));
    p_id = all_species_dict(char(protein_species(i)));
    dydt(mRNA_id) = transcription_rate(char(mRNA_species(i))) - (specific_growth_rate + transcription_correction(char(mRNA_species(i)))*degradation_constant_mRNA) * y(mRNA_id);
    dydt(p_id) = translation_rate(char(protein_species(i))) - (specific_growth_rate + translation_correction(char(protein_species(i)))*degradation_constant_protein) * y(p_id);
  end
  dydt = dydt.'
end