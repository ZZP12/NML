# specific to Julia, that is, related to the output file format (language).

#=
F: generate kinetic buffer
I:
O: return monod affinity constant symbol, mRNA species array,
  protein species array, W_string_array, disassociation_const_string_array.
BA
=#
function build_kinetics_buffer(all_species_dict::Dict{String, Int},
  all_rnx_list::Array, all_txtl_dict::Dict, sys2user::Dict)
  kinetics = "function calculate_kinetics(X, data_dictionary)"
  kinetics *= "\n\n\t# load all species dictionary"
  kinetics *= "\n\tall_species_dict = data_dictionary[\"all_species_dict\"]  # String-->Int"
  # for signaling
  kinetics *= "\n\n\t# generate kinetics for signaling, assume saturation kinetics"
  totalRnxNo = length(all_rnx_list)
  kinetics *= "\n\t# load Monod affinity constants"
  kinetics *= "\n\tMonodK = data_dictionary[\"MonodAffinityConstantDict\"]  # String-->Float"
  MonodAffinityConstant_String_Array = Array{String,1}()  # collect monod affinity constant symbol
  kinetics *= "\n\t# load reaction kinetic constants"
  kinetics *= "\n\tkcat = data_dictionary[\"kcat_signaling\"]  # in order of rnx"
  kinetics *= "\n\n\t# write reaction rate equations"
  kinetics *= "\n\trnx_rate_vector = zeros($totalRnxNo)"
  for (index, rnx) in enumerate(all_rnx_list)  # go thru every rnx
    tmp_rnx_rate = "kcat[$index]"
    if isdefined(rnx, :catalysts)  # enzyme
      for token in rnx.catalysts
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName])]"
      end
    end
    if isdefined(rnx,:reactants)  # substrate
      for token in rnx.reactants
        tmp_rnx_rate *= "*X[$(all_species_dict[token.oriBioName])]/(MonodK[\"MonodK_$(index)_$(token.oriBioName)\"]+X[$(all_species_dict[token.oriBioName])])"
        push!(MonodAffinityConstant_String_Array, "MonodK_$(index)_$(token.oriBioName)")
      end
    end
    kinetics *= "\n\trnx_rate_vector[$index] = "*tmp_rnx_rate
  end

  # for TXTL
  W_string_array = Array{String,1}()  # W: weight of protein action on transcription
  disassociation_const_string_array = Array{String,1}()  # disassociation constant in tranfer function
  kinetics *= "\n\n\t# generate kinetics & control terms for TXTL"
  kinetics *= "\n\t# load data from data dictionary"
  kinetics *= "\n\tW_value_dict = data_dictionary[\"W_value_dict\"]  #String-->Float"
  kinetics *= "\n\tcoop = data_dictionary[\"cooperativity\"]"
  kinetics *= "\n\tdisassociation_constant_dict = data_dictionary[\"transferFunctionDisassociationConstantDict\"]"
  kinetics *= "\n\tbackground_control_term = data_dictionary[\"backgroundGeneExpressionControlTermDict\"]"
  kinetics *= "\n\tkcatTX = data_dictionary[\"kcatTranscription\"]"
  kinetics *= "\n\tRNAPconc = data_dictionary[\"RNAPConcentration\"]"
  kinetics *= "\n\tgeneConc = data_dictionary[\"avgGeneConcentration\"]"
  kinetics *= "\n\tsaturationTX = data_dictionary[\"transcriptionSaturationConstant\"]"
  kinetics *= "\n\tkcatTL = data_dictionary[\"kcatTranslation\"]"
  kinetics *= "\n\tRIBOconc = data_dictionary[\"RIBOConcentration\"]"
  kinetics *= "\n\tsaturationTL = data_dictionary[\"translationSaturationConstant\"]"
  totalTXTLNo = length(all_txtl_dict)
  kinetics *= "\n\n\t# write control terms, transcription rate equations, and translation rate equations"
  kinetics *= "\n\tTX_control_term = Dict{String, Float64}()  # control term"
  kinetics *= "\n\tTX_rate_vector = Dict{String, Float64}()  # transcription rate"
  kinetics *= "\n\tTL_rate_vector = Dict{String, Float64}()  # translation rate"
  for (key, txtl) in all_txtl_dict # go thru every txtl
    tmp_protein_string = replace(key, sys2user["MRNA"], sys2user["PROTEIN"], 1)
    kinetics *= "\n\t# $key and $tmp_protein_string"
    up_factors_array = Array{String,1}()  # collection of upregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.activationProtein)  # upregulation
      for token_array in txtl.activationProtein
        tmp_up = ""
        tmp_up_name = key*"_"
        for token in token_array
          tmp_W = "W_$(key)_$(token.oriBioName)"
          tmp_up *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName])]^coop/((disassociation_constant_dict[\"KD_$(token.oriBioName)\"])^coop+X[$(all_species_dict[token.oriBioName])]^coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD_$(token.oriBioName)")
          tmp_up_name *= token.oriBioName*"_"
        end
        tmp_up = chop(tmp_up)
        tmp_up_name = chop(tmp_up_name)
        push!(up_factors_array, tmp_up_name)
        kinetics *= "\n\t$(tmp_up_name) = $(tmp_up)"
      end
    end
    down_factors_array = Array{String,1}()  # collection of downregulation factors name: targetedmRNA_factor(s)
    if !isempty(txtl.inhibitionProtein)  # downregulation
      for token_array in txtl.inhibitionProtein
        tmp_down = ""
        tmp_down_name = key*"_"
        for token in token_array
          tmp_W = "W_$(key)_$(token.oriBioName)"
          tmp_down *= "W_value_dict[\"$(tmp_W)\"] * X[$(all_species_dict[token.oriBioName])]^coop/((disassociation_constant_dict[\"KD_$(token.oriBioName)\"])^coop+X[$(all_species_dict[token.oriBioName])]^coop)*"
          push!(W_string_array, tmp_W)
          push!(disassociation_const_string_array, "KD_$(token.oriBioName)")
          tmp_down_name *= token.oriBioName*"_"
        end
        tmp_down = chop(tmp_down)
        tmp_down_name = chop(tmp_down_name)
        push!(down_factors_array, tmp_down_name)
        kinetics *= "\n\t$(tmp_down_name) = $(tmp_down)"
      end
    end
    # combine to form up and down, then merge up and down to form control term
    if !isempty(up_factors_array) && !isempty(down_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action))"
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action))"
      kinetics *= "\n\tTX_control_term[\"$key\"] = (background_control_term[\"$key\"] + up_action_on_$key)/(1 + background_control_term[\"$key\"] + up_action_on_$key + down_action_on_$key)"
    elseif !isempty(up_factors_array)
      up_action = ""
      for up in up_factors_array
        up_action *= up*" +"
      end
      kinetics *= "\n\tup_action_on_$key = $(chop(up_action))"
      kinetics *= "\n\tTX_control_term[\"$key\"] = (background_control_term[\"$key\"] + up_action_on_$key)/(1 + background_control_term[\"$key\"] + up_action_on_$key)"
    elseif !isempty(down_factors_array)
      down_action = ""
      for down in down_factors_array
        down_action *= down*" +"
      end
      kinetics *= "\n\tdown_action_on_$key = $(chop(down_action))"
      kinetics *= "\n\tTX_control_term[\"$key\"] = background_control_term[\"$key\"]/(1 + background_control_term[\"$key\"] + down_action_on_$key)"
    else
      # seems impossible, a TXTL dict is constructed b/c something acts on the key
    end
    kinetics *= "\n\tTX_rate_vector[\"$key\"] = TX_control_term[\"$key\"]*kcatTX*RNAPconc*geneConc/(saturationTX+geneConc)"
    kinetics *= "\n\tTL_rate_vector[\"$tmp_protein_string\"] = kcatTL*RIBOconc*X[$(all_species_dict[key])]/(saturationTL+X[$(all_species_dict[key])])"
  end

  kinetics *= "\n\n\treturn (rnx_rate_vector, TX_rate_vector, TL_rate_vector)"
  kinetics *= "\nend"
  return (kinetics, MonodAffinityConstant_String_Array, W_string_array, disassociation_const_string_array)
end

#=
F: read the file containing constants
I:
O:
BA
=#
function include_constants_from_literature(src_file_name, pad_string)
  # create src_buffer -
  src_buffer = ""
  # path to distrubtion -
  path_to_src_file = src_file_name
  open(path_to_src_file,"r") do src_file
    for line in eachline(src_file)
      src_buffer *= pad_string*"$line"
    end
  end
  return src_buffer
end

function build_data_dictionary_buffer(host_type::Symbol, all_species_array::Array,
  rnx_species_array::Array, all_rnx_list::Array, all_txtl_dict::Dict,
  Monod_affinity_constant_array::Array, W_string_array::Array,
  disassociation_const_array::Array, mRNA_species_array::Array, protein_species_array::Array)
  # data dictionary
  buffer = "function generate_model_parameters_dictionary()"

  # stoichiometry
  buffer *= "\n\n\t# load stoichiometry"
  buffer *= "\n\tstoichiometric_matrix = readdlm(\"./stoichiometry.dat\")"
  # all species dictionary
  buffer *= "\n\n\t# all species dictionary"
  buffer *= "\n\tall_species_dict = Dict{String, Int64}()  # String --> Int"
  #all_species_dict = Dict{String, Int64}()  # create all species dict: string --> int64
  for (index, st) in enumerate(all_species_array)
    buffer *= "\n\tall_species_dict[\"$st\"] = $index"
    #all_species_dict[st] = index
  end
  # all species dictionary_reversed
  buffer *= "\n\n\t# all species dictionary_reversed"
  buffer *= "\n\tall_species_reversed_dict = Dict{Int64, String}()  # Int --> String"
  buffer *= "\n\tfor (key,val) in all_species_dict"
  buffer *= "\n\t\tall_species_reversed_dict[val] = key"
  buffer *= "\n\tend"
  # all rnx species array
  buffer *= "\n\n\t# all rnx species array"
  buffer *= "\n\trnx_species_array::Array{String,1} = ["
  for st in rnx_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all mRNA species array
  buffer *= "\n\n\t# all mRNA species array"
  buffer *= "\n\tmRNA_species_array::Array{String,1} = ["
  for st in mRNA_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"
  # all protein species array
  buffer *= "\n\n\t# all protein species array"
  buffer *= "\n\tprotein_species_array::Array{String,1} = ["
  for st in protein_species_array
    buffer *= "\n\t\t\"$st\","
  end
  buffer = chop(buffer)*"\n\t]"

  # signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)
  buffer *= "\n\n\t# signaling reaction kinetic constants, reaction name: reactant(s)<rnx type:catalyst(s)>product(s)"
  buffer *= "\n\tkcat_signaling = 1e-3*ones($(length(all_rnx_list)))  # kcat[#reaction]: reaction name"
  for (id, rnx) in enumerate(all_rnx_list)
    buffer*= "\n\tkcat_signaling[$id] = 1.0  # kcat: $(rnx.rnxName)"
  end
  # Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate
  buffer *= "\n\n\t# Monod affinity constant in signaling reaction, nomenclature: MonodK_#reaction_Substrate"
  buffer *= "\n\tMonodAffinityConstantDict = Dict{String, Float64}()"
  for MK in Monod_affinity_constant_array
    # buffer *= "\n\t$MK = 1.0"
    # buffer *= "\n\tMonodAffinityConstant_dict[\"$MK\"] = $MK"
    buffer *= "\n\tMonodAffinityConstantDict[\"$MK\"] = 1.0"
  end

  # W_value in transcription control term, nomenclature: W_targetmRNA_actor
  buffer *= "\n\n\n\t# W_value in transcription control term, nomenclature: W_targetmRNA_actor"
  buffer *= "\n\tW_value_dict = Dict{String, Float64}()"
  for WS in W_string_array
    # buffer *= "\n\t$WS = 1.0"
    # buffer *= "\n\tW_value_dict[\"$WS\"] = $WS"
    buffer *= "\n\tW_value_dict[\"$WS\"] = 1.0"
  end
  # disassociation constants in transfer function, nomenclature: KD_speciesName
  buffer *= "\n\n\t# disassociation constants in transfer function, nomenclature: KD_speciesName"
  buffer *= "\n\ttransferFunctionDisassociationConstantDict = Dict{String, Float64}()"
  for DC in disassociation_const_array
    # buffer *= "\n\t$DC = 1.0"
    # buffer *= "\n\ttransferFunctionDisassociationConstantDict[\"$DC\"] = $DC"
    buffer *= "\n\ttransferFunctionDisassociationConstantDict[\"$DC\"] = 1.0"
  end
  # transcription specific correction factor
  buffer *= "\n\n\t# transcription specific correction factor"
  buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict = Dict{String, Float64}()"
  for mrna in mRNA_species_array
    buffer *= "\n\ttranscriptionSpecificCorrectionFactorDict[\"$mrna\"] = 1.0  # $mrna"
  end
  # background gene expression control term
  buffer *= "\n\n\t# background gene expression control term"
  buffer *= "\n\tbackgroundGeneExpressionControlTermDict = Dict{String, Float64}()"
  for mrna in mRNA_species_array
    buffer *= "\n\tbackgroundGeneExpressionControlTermDict[\"$mrna\"] = 0.1  # $mrna"
  end
  # translation specific correction factor
  buffer *= "\n\n\t# translation specific correction factor"
  buffer *= "\n\ttranslationSpecificCorrectionFactor = Dict{String, Float64}()"
  for pro in protein_species_array
    buffer *= "\n\ttranslationSpecificCorrectionFactor[\"$pro\"] = 1.0  # $pro"
  end
  # cooperativity number in transfer function
  buffer *= "\n\n\t# cooperativity number in transfer function"
  buffer *= "\n\tcooperativity = 1"
  #buffer *= "\n\n\tbackground_mRNA_synthesis_rate_vector = 0.01*ones($length_TXTL)"
  # load txtl constants buffer
  buffer *= "\n\n\n"
  if host_type == :bacteria
    buffer *= include_constants_from_literature("txtl_constants_ecoli.jl","\t")
  else
    buffer *= include_constants_from_literature("txtl_constants_hl60.jl","\t")
  end

  #---------------------------------
  # put all stuff in a Dictionary
  buffer *= "\n\n\t####################################"
  buffer *= "\n\t# put all stuff in a Dictionary"
  buffer *= "\n\tdataDictionary = Dict{String, Any}()"
  buffer *= "\n\tdataDictionary[\"stoichiometric_matrix\"] = stoichiometric_matrix"
  buffer *= "\n\tdataDictionary[\"all_species_dict\"] = all_species_dict"
  buffer *= "\n\tdataDictionary[\"all_species_reversed_dict\"] = all_species_reversed_dict"
  buffer *= "\n\tdataDictionary[\"rnx_species_array\"] = rnx_species_array"
  buffer *= "\n\tdataDictionary[\"mRNA_species_array\"] = mRNA_species_array"
  buffer *= "\n\tdataDictionary[\"protein_species_array\"] = protein_species_array"

  buffer *= "\n\tdataDictionary[\"kcat_signaling\"] = kcat_signaling"
  buffer *= "\n\tdataDictionary[\"MonodAffinityConstantDict\"] = MonodAffinityConstantDict"
  buffer *= "\n\tdataDictionary[\"W_value_dict\"] = W_value_dict"
  buffer *= "\n\tdataDictionary[\"transferFunctionDisassociationConstantDict\"] = transferFunctionDisassociationConstantDict"
  buffer *= "\n\tdataDictionary[\"transcriptionSpecificCorrectionFactorDict\"] = transcriptionSpecificCorrectionFactorDict"
  buffer *= "\n\tdataDictionary[\"backgroundGeneExpressionControlTermDict\"] = backgroundGeneExpressionControlTermDict"
  buffer *= "\n\tdataDictionary[\"translationSpecificCorrectionFactor\"] = translationSpecificCorrectionFactor"
  buffer *= "\n\tdataDictionary[\"cooperativity\"] = cooperativity"
  buffer *= "\n\tdataDictionary[\"RNAPConcentration\"] = rnapII_concentration"
  buffer *= "\n\tdataDictionary[\"RIBOConcentration\"] = ribosome_concentration"
  buffer *= "\n\tdataDictionary[\"degradationConstantmRNA\"] = degradation_constant_mRNA"
  buffer *= "\n\tdataDictionary[\"degradationConstantProtein\"] = degradation_constant_protein"
  buffer *= "\n\tdataDictionary[\"kcatTranscription\"] = kcat_transcription"
  buffer *= "\n\tdataDictionary[\"kcatTranslation\"] = kcat_translation"
  buffer *= "\n\tdataDictionary[\"avgGeneConcentration\"] = avg_gene_concentration"
  buffer *= "\n\tdataDictionary[\"transcriptionSaturationConstant\"] = saturation_transcription"
  buffer *= "\n\tdataDictionary[\"translationSaturationConstant\"] = saturation_translation"
  buffer *= "\n\tdataDictionary[\"specificGrowthRate\"] = maximum_specific_growth_rate - death_rate_constant"

  buffer *= "\n\n\treturn dataDictionary"
  buffer *= "\n\nend"

  return buffer
end


function build_simulation_buffer(NoExtracellularSpecies::Int64)
  buffer = build_copyright_header_buffer()
  buffer *= "\n\n# set up ODE & DAE"
  buffer *= "\nfunction simulation_ODEs(t,y,dataDictionary)"
  buffer *= "\n\t# get the derivatives"
  buffer *= "\n\n\t# correct for negatives -"
  buffer *= "\n\tidx_small = find(y.<0)"
  buffer *= "\n\ty[idx_small] = 0.0"

  buffer *= "\n\n\t# load data from data dictionary"
  buffer *= "\n\trnx_species::Array{String,1} = dataDictionary[\"rnx_species_array\"]"
  buffer *= "\n\tmRNA_species::Array{String,1} = dataDictionary[\"mRNA_species_array\"]"
  buffer *= "\n\tprotein_species::Array{String,1} = dataDictionary[\"protein_species_array\"]"
  buffer *= "\n\tstoichiometry = dataDictionary[\"stoichiometric_matrix\"]"
  buffer *= "\n\tall_species_dict::Dict{String, Int} = dataDictionary[\"all_species_dict\"]  # String --> number"
  buffer *= "\n\tspecific_growth_rate = dataDictionary[\"specificGrowthRate\"]"
  buffer *= "\n\ttranscription_correction = dataDictionary[\"transcriptionSpecificCorrectionFactorDict\"]"
  buffer *= "\n\ttranslation_correction = dataDictionary[\"translationSpecificCorrectionFactor\"]"
  buffer *= "\n\tdegradation_constant_mRNA = dataDictionary[\"degradationConstantmRNA\"]"
  buffer *= "\n\tdegradation_constant_protein = dataDictionary[\"degradationConstantProtein\"]"

  buffer *= "\n\n\t# kinetics"
  buffer *= "\n\t(rnx_rate, transcription_rate, translation_rate) = calculate_kinetics(y, dataDictionary)"

  buffer *= "\n\n\t# calculate the derivatives"
  buffer *= "\n\tdydt = zeros(length(all_species_dict))"
  buffer *= "\n\n\t# signaling reaction network"
  buffer *= "\n\tX2 = zeros(length(rnx_species)-$(NoExtracellularSpecies))  # species inside the cell"
  buffer *= "\n\tfor (id, st) in enumerate(rnx_species[$(NoExtracellularSpecies)+1:end])  # get X2 value from y"
  buffer *= "\n\t\tX2[id] = y[all_species_dict[st]]"
  buffer *= "\n\tend"
  buffer *= "\n\tdX1 = stoichiometry[1:$(NoExtracellularSpecies), :]*rnx_rate  # extracellular metabolites"
  buffer *= "\n\tdX2 = stoichiometry[$(NoExtracellularSpecies)+1:end, :]*rnx_rate - specific_growth_rate*X2"
  buffer *= "\n\tdX = [dX1; dX2]"
  buffer *= "\n\tfor (id, st) in enumerate(rnx_species)  # return derivatives back to dydt"
  buffer *= "\n\t\tdydt[all_species_dict[st]] = dX[id]"
  buffer *= "\n\tend"
  buffer *= "\n\n\t# TXTL network"
  buffer *= "\n\tfor i = 1:length(mRNA_species)  # seems easy to add additional terms to account for the signaling effect"
  buffer *= "\n\t\tmRNA_id = all_species_dict[mRNA_species[i]]"
  buffer *= "\n\t\tp_id = all_species_dict[protein_species[i]]"
  buffer *= "\n\t\tdydt[mRNA_id] = transcription_rate[mRNA_species[i]] - (specific_growth_rate + transcription_correction[mRNA_species[i]]*degradation_constant_mRNA)*y[mRNA_id]"
  buffer *= "\n\t\tdydt[p_id] = translation_rate[protein_species[i]] - (specific_growth_rate+ translation_correction[protein_species[i]]*degradation_constant_protein)*y[p_id]"
  buffer *= "\n\tend"
  buffer *= "\n\n\treturn dydt  # derivative"
  buffer *= "\nend"
  buffer *= "\n\n\nfunction simulation_DAEs(t, y, yp, r, dataDictionary)"
  buffer *= "\n\tdydt = simulation_ODEs(t, y, dataDictionary)"
  buffer *= "\n\tfor i = 1:length(dydt)"
  buffer *= "\n\t\tr[i] = dydt[i] - yp[i]"
  buffer *= "\n\tend"
  buffer *= "\nend"

  return buffer
end


function build_solveODEBalances_buffer(all_species_array::Array,
  all_species_dict::Dict{String, Int}, mRNA_species::Array, protein_species::Array)

  buffer = "include(\"./include.jl\")"

  buffer *= "\n\ntime_start = 0.0"
  buffer *= "\ntime_step_size = 0.1"
  buffer *= "\ntime_stop = 20.0"
  buffer *= "\nt = collect(time_start:time_step_size:time_stop)"

  No_species = length(all_species_dict)
  # initial concentration
  buffer *= "\n\n# initial species concentration"
  buffer *= "\ny0 = ["
  for i = 1:No_species-1
    buffer *= "\n\t0.0,  # $(all_species_array[i])"
  end
  buffer *= "\n\t0.0  # $(all_species_array[No_species])"
  buffer *= "\n]"

  # data dataDictionary
  buffer *= "\n\ndataDictionary = generate_model_parameters_dictionary()"
  buffer *= "\nall_species_reversed_dict = dataDictionary[\"all_species_reversed_dict\"]"
  # run simulation
  buffer *= "\n\nf(t, y) = simulation_ODEs(t, y, dataDictionary)"
  buffer *= "\nt, y = ode23s(f, y0, t)"
  # data transfer
  buffer *= "\n\n# data transfer"
  buffer *= "\nrow = length(t)"
  buffer *= "\ncol = length(y0)"
  buffer *= "\nY = zeros(row, col)"
  buffer *= "\nfor i = 1:col"
  buffer *= "\n\ttmp = map(x->x[i], y)"
  buffer *= "\n\tfor j = 1:row"
  buffer *= "\n\t\tY[j,i] = tmp[j]"
  buffer *= "\n\tend"
  buffer *= "\nend"
  # plotting
  buffer *= "\n\n# plot results"
  buffer *= "\nsubplt_col = round(Int, col/4)"
  buffer *= "\nfor i = 1:col"
  buffer *= "\n\tsubplot(4, subplt_col, i)"
  buffer *= "\n\tplot(t, Y[:,i])"
  buffer *= "\n\ttitle(all_species_reversed_dict[i])"
  buffer *= "\nend"

  return buffer
end

function build_copyright_header_buffer()
  current_year = string(Dates.year(now()))
  # Get comment data from
  buffer = ""
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  buffer*= "# Copyright (c) $(current_year) Varnerlab\n"
  buffer*= "# Robert Frederick Smith School of Chemical and Biomolecular Engineering\n"
  buffer*= "# Cornell University, Ithaca NY 14850\n"
  buffer*= "#\n"
  buffer*= "# Permission is hereby granted, free of charge, to any person obtaining a copy\n"
  buffer*= "# of this software and associated documentation files (the \"Software\"), to deal\n"
  buffer*= "# in the Software without restriction, including without limitation the rights\n"
  buffer*= "# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n"
  buffer*= "# copies of the Software, and to permit persons to whom the Software is\n"
  buffer*= "# furnished to do so, subject to the following conditions:\n"
  buffer*= "#\n"
  buffer*= "# The above copyright notice and this permission notice shall be included in\n"
  buffer*= "# all copies or substantial portions of the Software.\n"
  buffer*= "#\n"
  buffer*= "# THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n"
  buffer*= "# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n"
  buffer*= "# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n"
  buffer*= "# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n"
  buffer*= "# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n"
  buffer*= "# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN\n"
  buffer*= "# THE SOFTWARE.\n"
  buffer*= "# ----------------------------------------------------------------------------------- #\n"
  return buffer
end
