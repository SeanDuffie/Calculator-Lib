############################
# Imports
import math
#import matplotlib as mpl
#import matplotlib.pyplot as plt
############################


############################
# Define constants
SPD_LGT = 300000000 #in m/s
U_PERM = 0.00000047 #in H/m			Inductors
E_FRSPCE = 0.000000000008854 #in 	Capacitors
# 24x10 res
############################


############################
# Menu
def menu():
    """
    # Allows user to select through the options
    """
    a = ''
    while a != '.':
        print("Syllabus:")
        print("1. Real Components")
        print("2. Parasitic Effects")
        print("3. Pulse Waves")
        print("4. SS AC")
        print("5. Smith Chart")
        print("6. Microstrips")
        print("7. Waveguide")
        print("8. Antenna")
        print()

############################


############################
# Functions
"""

# Homework 1
def max_freq(lam):
    """
    #max_freq() -> float
    #Maximum frequency without propagation effects
    
    #Keyword arguments:
    #lam     -- length
"""
    f = SPD_LGT/(lam*10)
    print("Max Freq =", f, "Hz")
    return f

def inductance_coil(N, A, l):
    """
    # Inductance of a coil (Area is known)
"""
    L = (U_PERM*(N**2)*A)/l
    print("Inductance of Coil =", L, "H")
    return L

def capacitance_plates(A, d):
    """
    # Capacitance between two plates (Area is known)
"""
    C = (E_FRSPCE*A)/d
    print("Capacitance of Plates =", C, "F")
    return C

def res_freq(L, C):
    """
    # Resonance Frequency
"""
    fr = 1/(2*math.pi*math.sqrt(L*C))
    print("Resonance Frequency =", fr, "Hz")
    return fr

def ind_imp(f, L):
    """
    # Inductive Impedance
"""
    Z = (2*math.pi*f*L)
    print("Inductive Impedance =", Z, "Ohms")
    return Z

def par_ind(fr, C):
    """
    # Parasitic Inductance
"""
    L = ((1/(2*math.pi*fr))**2)/C
    print("Parasitic Inductance =", L, "H")
    return L

def cap_imp(f, C):
    """
    # Capacitive Impedance
"""
    Z = 1/(2*math.pi*f*C)
    print("Capacitive Impedance =", Z, "Ohms")
    return Z

def par_cap(fr, L):
    """
    # Parasitic Capacitance
"""
    C = ((1/(2*math.pi*fr))**2)/L
    print("Parasitic Capacitance =", C, "F")
    return C

def char_imp(L, C):
    """
    # Lossless Characteristic Impedance
"""
    Z = math.sqrt(L/C)
    print("Characteristic Impedance =", Z, "Ohms")
    return Z

def char_imp_l(L, R, C, G):
    """
    # Lossy Characteristic Impedance
"""
    Z = math.sqrt(L/C)
    print("Characteristic Impedance =", Z, "Ohms")
    return Z
"""
####################################################################################
# def plot_capimp_freq(C, f1, f2, s):	# Plots Capacitive Inductance over Frequency
#     X = math.linspace(f1, f2, s) # Parameters: (Start/Stop/Step)
#     Y = math.cap_imp(X, C)
    
#     fig, ac = plt.subplots()
#     ax.plot(X, Y, color='green')

#     fig.savefig("%.2fcap_%.0fto%.0fHz.pdf", C, f1, f2)
#     fig.show()

# def plot_indimp_freq(L, f1, f2, s):	# Plots Inductive Inductance over Frequency
#     X = math.linspace(f1, f2, s) # Parameters: (Start/Stop/Step)
#     Y = math.cap_imp(X, L)
    
#     fig, ac = plt.subplots()
#     ax.plot(X, Y, color='green')

#     fig.savefig("%.2find_%.0fto%.0fHz.pdf", L, f1, f2)
#     fig.show()
#####################################################################################
"""
def ref_coe(r0, r1):
    """
    # Finds the forward Reflection Coefficient
"""
    G = (r1-r0)/(r1+r0)
    print("Forward Reflection Coefficient =", G)
    return G
def ref_coe_l(r0, r2):
    """
    # Finds the reverse Reflection Coefficient
"""
    G = (r2-r0)/(r2+r0)
    print("Reverse Reflection Coefficient =", G)
    return G

def ref_coe_v1(v1, v2, r0):
    """
    # Uses V1, V2, and R0
"""
    print("Finding Reflection Coefficient 1...", v2/v1)
    return v2/v1
def ref_coe_v2(v1, v2, r0):
    """
    # Uses V1, V2, and R0
"""
    print("Finding Reflection Coefficients")
    return

#####################################################################################
def prop_delay_dv(d, v):
    """
    # Finds Prop Delay using VOP
"""
    t = (v/d)
    print("Prop Delay =", t, "seconds")
    return t

def prop_delay_dvf(d, vf):
    """
    # Finds Prop Delay using vf
"""
    v = vop_vf(vf)
    return prop_delay_dv(d, v)
def prop_delay_der(d, er):
    """
    # Finds Prop Delay using er
"""
    vf = 1/sqrt(er)
    print("Velocity Factor =", vf)
    return prop_delay_dvf(d, vf)
def vop_vf(vf):
    """
    # Finds VOP using VF
"""
    print("Velocity of Propagation =", SPD_LGT*vf, " m/s")
    return SPD_LGT*vf

def bounce_diagram1(vi, d, vf, r0, r1, r2):
    """
    # Constructs bounce diagrams using Volts in, cable distance, and VF
"""
    print("Constructing Bounce Diagrams, initial voltage")
    G1 = ref_coe(r0, r1)
    G2 = ref_coe_l(r0, r2)
    vr = vi*(r0/(r0 + r1))
    print("Initial Reflection =", vr, "V")
    while abs(vr) > 0.001:
        vr = vr*G2
        print("        Return Reflection =", vr, "V")
        vr = vr*G1
        print("    Forward Reflection =", vr, "V")
"""
############################
if __name__ == '__main__':
    menu()
    # max_freq(0.08)
    # ind = inductance_coil(16, 0.00003848, 0.012)
    # par_cap(82000000, ind)

menu()