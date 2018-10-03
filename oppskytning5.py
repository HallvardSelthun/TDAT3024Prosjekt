import numpy as np
from saturn_v import masse, skyvekraft, periode

# Konstanter for gravitasjon
radius_jord = 6371010
masse_jord = 5.9736 * 10**24
G_k = 6.67 * 10**(-11)


# Konstanter for luftmotstand
C_d = 0.25
areal_del_1 = np.pi * 5.55**2
areal_del_2 = np.pi * 5.55**2
areal_del_3 = np.pi * 3.3**2

print(areal_del_1, areal_del_2, areal_del_3)

# Tyngdekraft på raketten gitt av avstand til jorda og rakettens masse
def tyngdekraft(a, m):
    return G_k * masse_jord * m / (a**2)

# Atmosfærisk lag (1-4) som funksjon av høyde
def atm_lag(h):
    if(h<11000):
        return 1
    elif(h<25000):
        return 2
    elif(h<100000):
        return 3
    else:
        return 4

# Atmosfærens temperatur som funksjon av høyde
def temperatur(h):
    if(atm_lag(h)==1):
        return 288.19 - 0.00649*h
    elif(atm_lag(h)==2):
        return 216.69
    elif(atm_lag(h)==3):
        return 141.94 + 0.00299*h
    else:
        return 2.7

# Atmosfærens trykk som funksjon av høyde
def trykk(h):
    if(atm_lag(h) == 1):
        return 101290 * (temperatur(h)/288.08)**5.256
    elif(atm_lag(h) ==2):
        return 127760 * np.exp(-0.000157*h)
    elif(atm_lag(h) == 3):
        return 2488 * (temperatur(h)/216.6)**(-11.388)
    else:
        return 0


# Atmosfærens tetthet som funksjon av høyde
def tetthet(h):
    return 0.0034855*trykk(h)/temperatur(h)

# Luftmotstand som funksjon av rakettens høyde, hastighet og tid etter oppskytning
def luftmotstand(h0, v, t):
    h = (h0 - 6371000)
    # print("Luftmotstand: {}".format(0.5 * C_d * tetthet(h) * areal_del_1 * v**2))
    if(periode(t) == 0 or periode(t) == 1 or periode(t == 2)):
        return 0.5 * C_d * tetthet(h) * areal_del_1 * v**2
    else:
        return 0.5 * C_d * tetthet(h) * areal_del_3 * v**2



# for i in range(103):
#     print("temperatur: {}, trykk: {}, tetthet: {}".format(temperatur(i*1000), trykk(i*1000), tetthet(i*1000)))