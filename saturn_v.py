import numpy as np


# Konstanter(SI-enheter) for del 1
vekt_del_1 = 2290000
tom_vekt_del_1 = 131000
skyvekraft_del_1 = 35100000
brenntid_del_1 = 168

# Beregnede konstanter for del 1
vekttap_del_1 = vekt_del_1 - tom_vekt_del_1
masseendring_per_sek_del_1 = vekttap_del_1/brenntid_del_1
eksoshastighet_del_1 = skyvekraft_del_1/masseendring_per_sek_del_1

# Konstanter for del 2
vekt_del_2 = 496200
tom_vekt_del_2 = 40100
skyvekraft_del_2 = 5141000
brenntid_del_2 = 360

# Beregnede konstanter for del 2
vekttap_del_2 = vekt_del_2 - tom_vekt_del_2
masseendring_per_sek_del_2 = vekttap_del_2/brenntid_del_2
eksoshastighet_del_2 = skyvekraft_del_2/masseendring_per_sek_del_2

# Konstanter for del 3
vekt_del_3 = 123000
tom_vekt_del_3 = 13500
skyvekraft_del_3 = 1000000
brenntid_del_3 = 500

# Beregnede konstanter for del 3
vekttap_del_3 = vekt_del_3 - tom_vekt_del_3
masseendring_per_sek_del_3 = vekttap_del_3/brenntid_del_3
eksoshastighet_del_3 = skyvekraft_del_3/masseendring_per_sek_del_3

# Totale vektkonstanter
total_vekt_del_1 = vekt_del_1 + vekt_del_2 + vekt_del_3
total_vekt_del_2 = vekt_del_2 + vekt_del_3
total_vekt_del_3 = vekt_del_3

# Tidskonstanter
starttid_del_1 = 0
starttid_del_2 = brenntid_del_1
starttid_del_3 = brenntid_del_1 + brenntid_del_2
sluttid = brenntid_del_1 + brenntid_del_2 + brenntid_del_3


# Funksjon som returnerer 0-4 basert på hvilken forbrenningsperiode raketten befinner seg i
def periode(tid):
    if (tid<starttid_del_1):
        return 0
    elif(tid<=starttid_del_2):
        return 1
    elif(tid<=starttid_del_3):
        return 2
    elif(tid<=sluttid):
        return 3
    else:
        return 4

# Rakettens totale masse t sekunder etter oppskytning
def masse(tid):
    if(periode(tid)==0):
        return total_vekt_del_1
    elif(periode(tid)==1):
        return total_vekt_del_1 - masseendring_per_sek_del_1*(tid - starttid_del_1)
    elif(periode(tid)==2):
        return total_vekt_del_2 - masseendring_per_sek_del_2*(tid - starttid_del_2)
    elif(periode(tid)==3):
        return total_vekt_del_3 - masseendring_per_sek_del_3*(tid-starttid_del_3)
    else:
        return tom_vekt_del_3


def skyvekraft(tid):
    if (periode(tid)==0):
        return 0
    elif (periode(tid)==1):
        return skyvekraft_del_1
    elif (periode(tid)==2):
        return skyvekraft_del_2
    elif (periode(tid)==3):
        return skyvekraft_del_3
    else:
        return 0

