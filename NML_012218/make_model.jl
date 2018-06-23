using ArgParse
using Dates


function parse_commandline()
  settings_object = ArgParseSettings()
  @add_arg_table settings_object begin
    "-o"
      help = "Directory where the Julia model files will be written."
      arg_type = AbstractString
      default = "./autogeneratedmodel"

    "-m"
      help = "Path to the model file written in the NML format."
      arg_type = AbstractString
      required = true

    "-s"
      help = "Host type: bacteria or mammalian?"
      arg_type = Symbol
      default = :bacteria

    "-l"
      help = "Language: julia or python2 or python3 or matlab?"
      arg_type = AbstractString
      default = "julia"
  end
  return parse_args(settings_object)
end

parsed_args = parse_commandline()


include("preprocessor3.jl")
println("\n---------loading sentences------------")
inputSentencesArray = getRidOfNewlineHashcommentBlankline(parsed_args["m"])
tokenizedInputSentencesArray = sentenceTokenization2(inputSentencesArray)

println("\n------------normalizing and tagging------------")
reservedWordsPath = "reservedWords.jl"
taggedSentencesArray = tokenClassification2(tokenizedInputSentencesArray, reservedWordsPath)
# reshape for output observation
printTagSen = [["$(y[1])/$(y[2])" for y in x] for x in taggedSentencesArray]
println(typeof(printTagSen))
foreach(println, [join(ts, "  ") for ts in printTagSen])

include("biosymDecoder.jl")
println("\n------------information extraction------------")
include(reservedWordsPath)
BioSymVerbInfo = extractBioSymVerbInformation(taggedSentencesArray, reservedWords["SentenceType"])
# println("print by foreach:")
foreach(println, BioSymVerbInfo)

println("\n------------decoding bio symbols------------")
decodingBioSymGroups(BioSymVerbInfo)
println("\n------------decoding bio symbols results------------")
printArrayOfTripleToken(BioSymVerbInfo, AbstractString)
println("\n---------type conversion dictionary---------------")
typeConversionDict = setUpSymbolConvertionDict(BioSymVerbInfo)
foreach(println, typeConversionDict)
println("\n------replace each biosymbol string with generalBioSym composite-----------------")
newVerbBiosymInfo = replaceEachBioSymWithGeneralBioSym(BioSymVerbInfo, typeConversionDict)
printArrayOfTripleToken(newVerbBiosymInfo, generalBioSym)

println("\n------right before IR generation-----------------")
include("semanticChecking.jl")
preIR_DS = semanticCheckingForEachTriplet(newVerbBiosymInfo, typeConversionDict)
printArrayOfTripleToken(preIR_DS, generalBioSym)

println("\n-------model generation------------------")
include("modelGeneration.jl")
sys2userDict = Dict()
for (key, val) in typeConversionDict
  sys2userDict[val] = key
end
# reform to match with NML_V1
(rnx, txtl, rnx_set, mRNA_set, m_protein_set) = preparing_rnxList_txtlDict(
  preIR_DS, sys2userDict, typeConversionDict)
println("-----rnx list----")
foreach(x->println(x.rnxName), rnx)
# Sorting, insert "BIOMASS" here
(sorted_rnx_species_array, rnx_species2index_dict, rnx_index2species_dict,
 sorted_all_species_array, all_species2index_dict, all_index2species_dict, extra_species_num) =
  sorting_species_list(rnx_set, mRNA_set, m_protein_set, sys2userDict)
println("------------")
foreach(println, sorted_rnx_species_array)
println("------------")
foreach(println, sorted_all_species_array)

# write program to disk
output_file = Dict{String, String}()  # put all files in this dict
# stoichiometric_matrix
stoichiometric_matrix_buffer = build_stoichiometric_matrix_buffer(rnx, rnx_species2index_dict)
output_file["stoichiometry.dat"] = stoichiometric_matrix_buffer


target_lang = parsed_args["l"]
if target_lang == "julia"
  include("JuliaStrategy.jl")
  file_suffix = ".jl"
elseif target_lang == "python3"
  include("PythonStrategy.jl")
  file_suffix = ".py"
elseif  target_lang == "matlab"
  include("MATLABStrategy.jl")
  file_suffix = ".m"
else
  println("Generate model in julia as default")
  include("JuliaStrategy.jl")
  file_suffix = ".jl"
end

(kinetics_buffer, Monod_const, W_array, disassociation_const) =
  build_kinetics_buffer(all_species2index_dict, rnx, txtl, sys2userDict)
output_file["calculate_kinetics" * file_suffix] = kinetics_buffer

host_type = parsed_args["s"]
data_dictionary_buffer = build_data_dictionary_buffer(host_type, sorted_all_species_array,
  all_species2index_dict,
  sorted_rnx_species_array, rnx, txtl, Monod_const, W_array, disassociation_const,
  collect(mRNA_set), collect(m_protein_set))
output_file["generate_dataDictionary" * file_suffix] = data_dictionary_buffer

ODE_simulation_buffer = build_simulation_buffer(extra_species_num)
output_file["ODEbalance" * file_suffix] = ODE_simulation_buffer

solveODEBalances_buffer = build_solveODEBalances_buffer(sorted_all_species_array,
  all_species2index_dict, collect(mRNA_set), collect(m_protein_set))
output_file["solveODEbalance" * file_suffix] = solveODEBalances_buffer


current_dir = parsed_args["o"]
if !(isdir(current_dir))
  mkpath(current_dir)
end
if target_lang == "julia"
  cp("include.jl", current_dir*"/include.jl"; remove_destination=true)
end
for (file_name, file) in output_file
  write(current_dir*"/"*file_name, file)
end