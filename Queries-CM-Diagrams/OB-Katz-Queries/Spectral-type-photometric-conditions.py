import numpy as np
import sys
import os
import imageio.core
from astropy.io import fits

spectraltypeO = np.genfromtxt("Data/SpectraltypeO.dat", names = True, dtype=None)
VminIc_O = spectraltypeO['VIc']
V_O = spectraltypeO['Mv']
HminKs_O = spectraltypeO['HKs']

G_O = V_O - 0.01746 + 0.008092*(VminIc_O) - 0.2810*(VminIc_O)**2 + 0.03655*(VminIc_O)**3
G_bpminG_rp_O = -0.04212 + 1.286*(VminIc_O) - 0.09494*(VminIc_O)**2
G_bpminG_rp_O_2mass = 0.4238 + 6.098*(HminKs_O) + 1.991*(HminKs_O)**2

print("Conditions for Spectral type O are:")
print("-infinity < M_G < "+str(np.max(G_O)))
print("Johnson-Cousins relationship: -infinity < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_O)))
print("2MASS relationship: -infinity < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_O_2mass)))
print("")

spectraltypeB = np.genfromtxt("Data/SpectraltypeB.dat", names = True, dtype=None)
VminIc_B = spectraltypeB['VIc']
V_B = spectraltypeB['Mv']
HminKs_B = spectraltypeB['HKs']

G_B = V_B - 0.01746 + 0.008092*(VminIc_B) - 0.2810*(VminIc_B)**2 + 0.03655*(VminIc_B)**3
G_bpminG_rp_B = -0.04212 + 1.286*(VminIc_B) - 0.09494*(VminIc_B)**2
G_bpminG_rp_B_2mass = 0.4238 + 6.098*(HminKs_B) + 1.991*(HminKs_B)**2

print("Conditions for Spectral type B are:")
print(str(np.min(G_B))+" < M_G < "+str(np.max(G_B)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_B))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_B)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_B_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_B_2mass)))
print("")

spectraltypeA = np.genfromtxt("Data/SpectraltypeA.dat", names = True, dtype=None)
VminIc_A = spectraltypeA['VIc']
V_A = spectraltypeA['Mv']
HminKs_A = spectraltypeA['HKs']

G_A = V_A - 0.01746 + 0.008092*(VminIc_A) - 0.2810*(VminIc_A)**2 + 0.03655*(VminIc_A)**3
G_bpminG_rp_A = -0.04212 + 1.286*(VminIc_A) - 0.09494*(VminIc_A)**2
G_bpminG_rp_A_2mass = 0.4238 + 6.098*(HminKs_A) + 1.991*(HminKs_A)**2

print("Conditions for Spectral type A are:")
print(str(np.min(G_A))+" < M_G < "+str(np.max(G_A)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_A))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_A)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_A_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_A_2mass)))
print("")

spectraltypeF = np.genfromtxt("Data/SpectraltypeF.dat", names = True, dtype=None)
VminIc_F = spectraltypeF['VIc']
V_F = spectraltypeF['Mv']
HminKs_F = spectraltypeF['HKs']

G_F = V_F - 0.01746 + 0.008092*(VminIc_F) - 0.2810*(VminIc_F)**2 + 0.03655*(VminIc_F)**3
G_bpminG_rp_F = -0.04212 + 1.286*(VminIc_F) - 0.09494*(VminIc_F)**2
G_bpminG_rp_F_2mass = 0.4238 + 6.098*(HminKs_F) + 1.991*(HminKs_F)**2

print("Conditions for Spectral type F are:")
print(str(np.min(G_F))+" < M_G < "+str(np.max(G_F)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_F))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_F)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_F_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_F_2mass)))
print("")

spectraltypeG = np.genfromtxt("Data/SpectraltypeG.dat", names = True, dtype=None)
VminIc_G = spectraltypeG['VIc']
V_G = spectraltypeG['Mv']
HminKs_G = spectraltypeG['HKs']

G_G = V_G - 0.01746 + 0.008092*(VminIc_G) - 0.2810*(VminIc_G)**2 + 0.03655*(VminIc_G)**3
G_bpminG_rp_G = -0.04212 + 1.286*(VminIc_G) - 0.09494*(VminIc_G)**2
G_bpminG_rp_G_2mass = 0.4238 + 6.098*(HminKs_G) + 1.991*(HminKs_G)**2

print("Conditions for Spectral type G are:")
print(str(np.min(G_G))+" < M_G < "+str(np.max(G_G)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_G))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_G)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_G_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_G_2mass)))
print("")

spectraltypeK = np.genfromtxt("Data/SpectraltypeK.dat", names = True, dtype=None)
VminIc_K = spectraltypeK['VIc']
V_K = spectraltypeK['Mv']
HminKs_K = spectraltypeK['HKs']

G_K = V_K - 0.01746 + 0.008092*(VminIc_K) - 0.2810*(VminIc_K)**2 + 0.03655*(VminIc_K)**3
G_bpminG_rp_K = -0.04212 + 1.286*(VminIc_K) - 0.09494*(VminIc_K)**2
G_bpminG_rp_K_2mass = 0.4238 + 6.098*(HminKs_K) + 1.991*(HminKs_K)**2

print("Conditions for Spectral type K are:")
print(str(np.min(G_K))+" < M_G < "+str(np.max(G_K)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_K))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_K)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_K_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_K_2mass)))
print("")

spectraltypeM = np.genfromtxt("Data/SpectraltypeM.dat", names = True, dtype=None)
VminIc_M = spectraltypeM['VIc']
V_M = spectraltypeM['Mv']
HminKs_M = spectraltypeM['HKs']

G_M = V_M - 0.01746 + 0.008092*(VminIc_M) - 0.2810*(VminIc_M)**2 + 0.03655*(VminIc_M)**3
G_bpminG_rp_M = -0.04212 + 1.286*(VminIc_M) - 0.09494*(VminIc_M)**2
G_bpminG_rp_M_2mass = 0.4238 + 6.098*(HminKs_M) + 1.991*(HminKs_M)**2

print("Conditions for Spectral type M are:")
print(str(np.min(G_M))+" < M_G < "+str(np.max(G_M)))
print("Johnson-Cousins relationship: " + str(np.min(G_bpminG_rp_M))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_M)))
print("2MASS relationship: " + str(np.min(G_bpminG_rp_M_2mass))+" < G_{BP} - G_{RP} < "+str(np.max(G_bpminG_rp_M_2mass)))
