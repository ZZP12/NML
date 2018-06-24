# washout simulation
include("include.jl")
include("SolveBalances.jl")

# set up simulation time for each phase
time_start = 0
PI_dura = 7
PII_dura = 10
PIII_dura = 8
time_step = 0.1

# phase I: run to steady state
DataDict = generate_model_parameters_dictionary()
TIEnd = time_start + PI_dura
(TP1, XP1) = SolveBalances(time_start, TIEnd, time_step, DataDict)
# phase II: add substrate
IniCond = XP1[end, :]
IniCond[1] = 0.001
IniCond[2] = 10
DataDict["initial_condition"] = IniCond
TIIEnd = TIEnd + PII_dura
(TP2, XP2) = SolveBalances(TIEnd, TIIEnd, time_step, DataDict)
# phase III: remove substrate
IniCond = XP2[end, :]
IniCond[2] = 0
DataDict["initial_condition"] = IniCond
TIIIEnd = TIIEnd + PIII_dura
(TP3, XP3) = SolveBalances(TIIEnd, TIIIEnd, time_step, DataDict)

# pack the three phases together
T = [TP1; TP2; TP3]
X = [XP1 ; XP2 ; XP3];
# plot results
all_species_reversed_dict = DataDict["all_species_reversed_dict"]
subplt_col = ceil(Int, length(IniCond)/4)
plt.figure(tight_layout=true)
for i = 1:length(IniCond)
  plt.subplot(4, subplt_col, i)
  plt.plot(T[40:end]-4, X[40:end,i])
  plt.title(all_species_reversed_dict[i])
end
plt.show()

# writedlm("simulation.dat", [T X])

# column; species
# 1  BIOMASS
# 2  met_A_e
# 3  met_B_e
# 4  met_C_e
# 5  met_A_c
# 6  met_B_c
# 7  met_C_c
# 8  p_E1_c
# 9  p_E2_c
# 10  p_E3_c
# 11  p_E4_c
# 12  p_TA_c
# 13  p_TB_c
# 14  p_TC_c
# 15  mRNA_E1_c
# 16  mRNA_E2_c
# 17  mRNA_E3_c
# 18  mRNA_E4_c
# 19  mRNA_TA_c
# 20  mRNA_TB_c
# 21  mRNA_TC_c
