#
   #

# P_, p_, m_, g_, M_ are PROTEIN, PROTEIN, MRNA, GENE, METABOLITE
# _e, _c, -p are _exc, _c, -pho
# p_C & p_H reacting into [] irreversibly
# p_C irreversibly reacting into []
# (P_A or P_B and P_C promotes the expression of P_D
# (P_A or P_B) promotes the expression of P_D
# (P_A and P_B) promotes the expression of P_D

# the cell uptakes M_L_e
# the cell uptakes m_L via p_H   # does not check compartment
# the cell uptakes 2*m_L_e via p_H_c
# secretes m_L via p_L   # simple form
# the cell secretes m_J_c
# the cell secretes M_J_c
# 2*p_C and p_m binds together
# 3*p_C binding together into p_L
# 3*p_C  binds into p_L
#
# p_L unbinds
# p_L unbinding into 3*p_m
# p_A and p_B phosphorylate p_D or p_F
#
# p_A phosphorylating p_D at site S50
# p_A phosphorylating p_D
# p_A dephosphorylated p_F-p
# p_A dephosphorylated p_F-S50 at site S50
# p_A dephosphorylated p_H at site S60
#
# p_C activates the expression of g_B
# m_A inhibited the transcription of g_F
# p_L activated the translation of g_m
# p_D promote p_E
# m_B inhibiting g_C

# P_A, P_B, and P_C promotes P_C or P_D or P_F  # ambiguity
# (P_A or P_B) and P_C and p_D and (P_E or p_F) promotes the expression of g_E and g_F
# (P_A and P_B) or P_C or P_D or (P_E and P_F) promotes the expression of (P_A or P_B) and (P_E or p_F)  # ambiguity
# (P_A and P_B) or (P_E and p_F) promotes (P_A or P_B) and P_C  # ambiguity
# # sentences for testing
# (P_A and P_B) or P_C or P_D and (P_E and P_F) promotes the expression of g_E
# (P_A and P_B) or P_C or P_D or (P_E or P_F) promotes the expression of g_E

# (P_A or P_B) and P_C and p_D and (P_E or p_F) inhibit the expression of p_E and p_F

# p_F catalyzing the conversion of 3*P_F and 2*p_H into 2*p_H
# p_F catalyzing the reversible conversion of 3*P_F and 2*p_H into 2*p_H
# p_F catalyzing the conversion reversible of 3*P_F and 2*p_H into 2*p_H  # cool or funny
# p_F or p_C catalyzing 3*P_F and 2*p_H into 2*p_H
#
# p_F catalyzes 3*P_F and 2*p_H into 2*p_H
#
 # p_A pmomotes  the expression of p_C
