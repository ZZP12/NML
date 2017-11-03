#
   #

P_, p_, m_, g_, M_ are PROTEIN, PROTEIN, MRNA, GENE, METABOLITE
p_C reacting into p_D
# (P_A or P_B and P_C promotes the expression of P_D
# (P_A or P_B) promotes the expression of P_D
# (P_A and P_B) promotes the expression of P_D

# the cell uptakes m_L
# the cell uptakes m_L via p_H
# secretes m_L via p_L   # simple form
# the cell secretes m_J
# 2*p_C binds together
# 3*p_C binding together into p_L
# 3*p_C  binds into p_L
#
# p_L unbinds
# p_L unbinding into 3*p_m
# p_A and p_B phosphorylate p_D or p_F
#
# p_A phosphorylating p_D at site S50
# p_A dephosphorylated p_F
# p_A dephosphorylated p_H at site S60
#
# p_C activates the expression of g_B
# m_A inhibited the transcription of g_F
# p_L activated the translation of g_m
# p_D promote p_E
# m_B inhibiting g_C

# P_A, P_B, and P_C promotes P_C or P_D or P_F
# (P_A or P_B) and P_C and p_D and (P_E or p_F) promotes the expression of g_E
# (P_A and P_B) or P_C or P_D or (P_E and P_F) promotes the expression of (P_A or P_B) and (P_E or p_F)
# (P_A and P_B) or (P_E and p_F) promotes (P_A or P_B) and P_C
# # sentences for testing
# (P_A and P_B) or P_C or P_D and (P_E and P_F) promotes the expression of g_E
# (P_A and P_B) or P_C or P_D or (P_E or P_F) promotes the expression of g_E


# p_F catalyzing the conversion of 3*P_F and 2*p_H into 2*p_H
p_F catalyzing the reversible conversion of 3*P_F and 2*p_H into 2*p_H
# p_F catalyzing the conversion reversible of 3*P_F and 2*p_H into 2*p_H  # cool or funny
# p_F or p_C catalyzing 3*P_F and 2*p_H into 2*p_H
#
# p_F catalyzes 3*P_F and 2*p_H into 2*p_H
#
 # p_A pmomotes  the expression of p_C



# # simplified prostate cancer network
# # symbol conversion
# g_, m_, p_, met_, P_ are TAGgene, TAGmRNA, TAGprotein, TAGmetabolite, TAGprotein
# ph, e,  c are PHOS, EXTRA,  CYT
#
# # signaling
# The cell uptake met_DHEA_e & met_Testosterone_e
# [] catalyzes the conversion of met_DHEA_c into met_AD_c
# [] catalyzes the conversion of met_Testosterone_c into met_DHT_c
# [] catalyzes the conversion of met_AD_c into 0.1*met_Testosterone_c & 0.9*met_DHT_c
# [] catalyzes the conversion of 2*met_DHT_c & 2*p_AR_c into p_DHT_AR_ph_2_c
#
# [] catalyzes the reversible conversion of 2*p_Her2_c into p_Her2_ph_2_c
# p_Her2_ph_2_c phosphorylates p_PI3K_c
# p_PI3K_ph_c phosphorylates p_Akt_c
# p_PTEN_c dephosphorylates p_Akt_ph_c
# p_Her2_ph_2_c catalyzes the conversion of p_MAPK_c into p_MAPK_act_c
# p_MAPK_act_c catalyze the reversible conversion of 2*p_AR_c into p_AR_ph_2_c
# p_Akt_ph_c catalyze the conversion of 2*p_AR_c into p_AR_ph_2_c
#
# # TXTL
# p_MAPK_act_c activates p_AP1_c
# p_DHT_AR_ph_2_c & p_AR_ph_2_c activate p_PSA_c
# p_AR_ph_2_c activates p_AR_c
# p_AP1_c downregulates p_AR_c
#
# # degradation
# p_DHT_AR_ph_2_c, p_AR_c, p_PSA_c, p_AR_ph_2_c, p_Her2_c, p_PI3K_c, p_PTEN_c, p_Akt_c, p_MAPK_c, p_AP1_c react into []
