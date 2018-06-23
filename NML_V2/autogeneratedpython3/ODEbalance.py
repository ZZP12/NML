
# set up ODE, get the derivatives
def simulation_ODEs(y, t, dataDictionary):
  # from dataDictionary import generate_model_parameters_dictionary
  from kinetics import calculate_kinetics
  import numpy as np
  # dataDictionary = generate_model_parameters_dictionary()
  # correct for negatives
  y = [x if x > 0 else 0 for x in y]

  # load data from data dictionary
  rnx_species = dataDictionary["rnx_species_array"]
  mRNA_species = dataDictionary["mRNA_species_array"]
  protein_species = dataDictionary["protein_species_array"]
  stoichiometry = dataDictionary["stoichiometric_matrix"]
  all_species_dict = dataDictionary["all_species_dict"]  # String --> number
  specific_growth_rate = dataDictionary["specificGrowthRate"]
  transcription_correction = dataDictionary["transcriptionSpecificCorrectionFactorDict"]
  translation_correction = dataDictionary["translationSpecificCorrectionFactor"]
  degradation_constant_mRNA = dataDictionary["degradationConstantmRNA"]
  degradation_constant_protein = dataDictionary["degradationConstantProtein"]

  # kinetics
  rnx_rate, transcription_rate, translation_rate = calculate_kinetics(y, dataDictionary)

  # calculate the derivatives
  dydt = np.zeros(len(all_species_dict))

  # signaling reaction network
  X2 = np.zeros(len(rnx_species)-2)  # species inside the cell
  for id, st in enumerate(rnx_species[2:]):
    # get X2 value from y
    X2[id] = y[all_species_dict[st]]

  dX1 = np.dot(stoichiometry[0:2, :], rnx_rate)  # extracellular metabolites
  dX2 = np.dot(stoichiometry[2: , :], rnx_rate) - specific_growth_rate * X2  # INTRA
  # print(np.shape(dX1), np.shape(dX2))
  dX = np.concatenate((dX1, dX2), axis=1)
  # print(np.shape(dX))
  for id, st in enumerate(rnx_species):
    # return derivatives back to dydt
    dydt[all_species_dict[st]] = dX[0, id]

  dydt[0] = specific_growth_rate * y[0]  # for Biomass

  # TXTL network
  for i in range(len(mRNA_species)):
  # seems easy to add additional terms to account for the signaling effect
    mRNA_id = all_species_dict[mRNA_species[i]]
    p_id = all_species_dict[protein_species[i]]
    dydt[mRNA_id] = transcription_rate[mRNA_species[i]] - (specific_growth_rate + transcription_correction[mRNA_species[i]]*degradation_constant_mRNA) * y[mRNA_id]
    dydt[p_id] = translation_rate[protein_species[i]] - (specific_growth_rate + translation_correction[protein_species[i]]*degradation_constant_protein) * y[p_id]


  return dydt