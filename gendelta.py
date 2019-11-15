import numpy as np
import sys

n_spec = 0.96
ACDM = 7.79627709955e-09
z_eq = 3379.0
z_output = 8000000.0
Omega = 0.3
OmegaLambda = 0.7
HubbleParam = 0.7
UnitLength_in_cm = 3.085678e24  # Mpc
UnitVelocity_in_cm_per_s =  1e5 # km/s
kmin = 0.001
kmax = 1e10
def main():
    p = np.loadtxt(sys.argv[1])
    k = p[:,0]
    T = p[:,1]
    Delta2 = k**(n_spec+3)*T**2*ACDM/(2.*np.pi**2)
    #k_ref = k[len(k)-1]
    #Delta_ref = create_delta(k_ref,z_output)
    #Ratio = Delta2[len(k)-1]/Delta_ref

    #k_new = np.arange(np.log10(kmin),np.log10(kmax),0.001)
    #k_new = 10.0**k_new
    #Delta_new = Ratio*create_delta(k_new, z_output)    
    for i in range(len(k)):
        print np.log10(k[i]),np.log10(Delta2[i])

        
# def create_delta(kp,z):
#     UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s
#     Hubble = 3.2407789e-18 * UnitTime_in_s
#     OmegaRadiation = Omega/(2+z_eq)
#     a_eq = 1.0/(1+z_eq)
#     a = 1./(1+z)
#     y = a/a_eq
#     T_cmb = 2.7255 # // Fixsen 2009
#     N_eff = 3.046 # // Mangano et al (2002,2005)
#     T_ref = 2.7 #;
#     Theta = T_cmb/T_ref #;
#     f_v = 1.0 - 1.0/(N_eff*(7./8.)*(4./11.)**(4./3.) + 1.0) #;
#     I_2 = 0.594*(1 - 0.631*f_v + 0.284*f_v*f_v) #; // eq B14 Hu & Sukiyama 1996
#     k_eq  = 9.67e-2 * Omega * HubbleParam * HubbleParam * np.sqrt(1.0 - f_v) / (Theta*Theta) #; // 1/Mpc

#     k = kp * HubbleParam #; //change to 1/Mpc
#     hubble_a = (Hubble * np.sqrt(Omega / a**3 + OmegaRadiation/a**4 + 
#                               (1 - Omega - OmegaLambda ) / a**2 + OmegaLambda) )
#     y = a/a_eq #;
    
#     aHaEQ = (1.0 + np.sqrt(1.0 + 8.0*(k/k_eq)*(k/k_eq)))/(4.0 * (k/k_eq)*(k/k_eq)) #; //eq B13 Hu & Sukiyama 1996
#     logk_term = np.log(4.0 * I_2 * np.exp(-3.0) / aHaEQ) #;
#     d = (logk_term - np.log(( np.sqrt(1 + y) + 1 )/( np.sqrt(1 + y) - 1)))*(y + 2.0/3.0) + 2.0*np.sqrt(1 + y) #;  //eq D3 Hu & Sukiyama 1996
#     return d*d

if __name__ == "__main__":
    main()
