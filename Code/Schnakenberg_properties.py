# =================================================================================
# =================================================================================
# Script:"Schnakenberg_properties"
# Date: 2021-03-08
# Implemented by: Johannes Borgqvist
# Description:
# This script calculates the properties of the Schnakenberg model. 
# =================================================================================
# =================================================================================
# Import Libraries
# =================================================================================
# =================================================================================
# Import numpy as well
import numpy as np
# Write less trick
np.fac = np.math.factorial
# =================================================================================
# =================================================================================
# Functions for calculating the properties of the Schnakenberg model.
# =================================================================================
# =================================================================================

#------------------------------------------------------------------
# Function 1: "calculate_critical_parameters_Schnakenberg"
# The function takes the two parameters a and b of the
# Schnakenberg RD modelas well as the squared wavenumber
# k_squared as input. It returns the criticald diffusion
# parameter dc and the critical wavenumber gamma_c.
#------------------------------------------------------------------

def calculate_steady_states_and_critical_parameters_Schnakenberg(a,b,k_squared):

    # Calculate the steady states
    u_0 = a + b
    v_0 = ( (b) / ( (a + b)** 2))

    # Calculate the four partial derivatives of
    # the Jacobian matrix used in the linear
    # stability analysis
    f_u = ( (b - a) / (b + a) )
    f_v = ( (a + b)**2 )
    g_u = ( (-2*b) / (a+b) )
    g_v = -f_v

    # Calculate the critical diffusion value reported
    # in Chaplain 2001. Since the expression is quite
    # messy, we divide it into three parts.
    d_c = ( (f_u*g_v) - (2*f_v*g_u) )
    d_c +=  np.sqrt( d_c**2 - ((f_u**2)*(g_v**2)) )
    d_c = ( (d_c) / (f_u**2) )

    # Using the critical value of the diffusion d_c,
    # we also calculate the critical wave length number
    # gamma_c reported in Chaplain 2001.
    gamma_c = ( ( 2 * d_c * k_squared) / ( (d_c*f_u) + g_v ) )
    
    # Return the parameters 
    return u_0, v_0, d_c, gamma_c

#------------------------------------------------------------------
# Function 2: "compute_minimal_holeradius_for_pattern_disturbance"
# The function takes the Schnakenberg model parameters a, b, d, and
# gamma, as well as the wavenumber k_squared. It computes and
# returns the minimal radius epsilon of a geodesic disk-shaped hole
# that is needed to affect the Turing patterns. The function also
# returns the spectral parameters n_eps and m_eps that corresponds
# to the eigenfunction(s) causing the disturbance
#------------------------------------------------------------------

def compute_minimal_holeradius_for_pattern_disturbance(a,b,d,gamma,n_ref,n_largest):

    # Compute the lower and upper excited wavenumber bounds
    # OBS! There should be only one eigenvalue in the excitation
    # interval (gL, gM), namely n_ref*(n_ref+1)
    lm = d*(b-a) - (a+b)**3
    gL = gamma*(lm - np.sqrt(lm**2 -4*d*(a+b)**4))/(2*d*(a+b))
    gM = gamma*(lm + np.sqrt(lm**2 -4*d*(a+b)**4))/(2*d*(a+b))

    # Initialize list for epsilon, n, m
    eps_list = [] 
    
    # Go through LOWER unperturbed eigenvalue parameters to compute
    # the minimal epsilon needed to enter the excitation interval
    # For LOWER n's we enter from below, i.e., hit gL
    for n in range(1,n_ref):

        eigv_n = n*(n+1)          # The unperturbed eigenvalue
        m = 0                     # Corresponding and relevant m-value
        Cnm = 4*(2*n+1)/(n*(n+1)) # Coefficient of leading perturbation term     

        epsilon = np.sqrt((gL - eigv_n)/Cnm)
        eps_list.append([n, m, epsilon])
        #print(n, m, epsilon)


        
    # Go through the given unperturbed eigenvalue parameter n_ref to compute
    # the minimal epsilon needed to exit the excitation interval
    # For n_ref we can exit both below and above, i.e., hit gL and gM
    #print("----")

    n = n_ref  
    eigv_n = n*(n+1)          # The unperturbed eigenvalue

    # m = 0, we can only exit above, i.e., hit gM
    m = 0                     # Corresponding and relevant m-value
    Cnm = 4*(2*n+1)/(n*(n+1)) # Coefficient of leading perturbation term       

    epsilon = np.sqrt((gM - eigv_n)/Cnm)
    eps_list.append([n, m, epsilon])
    #print(n, m, epsilon)

    # m > 0, we can only exit below, i.e., hit gL
    # Loop over the corresponding and relevant m-values
    for m in range(1,n+1):

            Cnm_nom = (2*n+1)*np.fac(m+n)
            Cnm_den = (4**m)*np.fac(n-m)*np.fac(m)*np.fac(m-1)
            Cnm = -Cnm_nom/Cnm_den # Coefficient of leading perturbation term

            epsilon = ((gL - eigv_n)/Cnm)**(1/(2*m))
            eps_list.append([n, m, epsilon])
            #print(n, m, epsilon)
            
    #print("----")


    
    # Go through some HIGHER unperturbed eigenvalue parameters to compute
    # the minimal epsilon needed to enter the excitation interval
    # For HIGHER n's we enter from above, i.e., hit gM
    for n in range(n_ref+1,n_largest+1):
        eigv_n = n*(n+1)          # The unperturbed eigenvalue      
        # Loop over the corresponding and relevant m-values
        # (m = 1, ..., n for HIGHER n's)
        for m in range(1,n+1):
            Cnm_nom = (2*n+1)*np.fac(m+n)
            Cnm_den = (4**m)*np.fac(n-m)*np.fac(m)*np.fac(m-1)
            Cnm = -Cnm_nom/Cnm_den # Coefficient of leading perturbation term
            epsilon = ((gM - eigv_n)/Cnm)**(1/(2*m))
            eps_list.append([n, m, epsilon])
            #print(n, m, epsilon)



    # Compute the minimal epsilon and corresponding spectral parameters        
    #print(eps_list)
    eps_min = 3.14 # Initialize with something large.
    for i in range(0, len(eps_list)):
        #print(eps_list[i])
        if eps_list[i][2] < eps_min:
            eps_min = eps_list[i][2]
            n_min = eps_list[i][0]
            m_min = eps_list[i][1]

    #print(eps_min, n_min, m_min)
    
    # Return the parameters 
    return eps_min, n_min, m_min
#------------------------------------------------------------------
# Function 3: "perturbed_eigenvalues_Schnakenberg"
# This function calculates the perturbed eigenvalues of the
# Schnackenberg model. It takes the eigen value pair (n,m)  as well
# as the hole radius epsilon and return the eigenvalue. 
#------------------------------------------------------------------
def perturbed_eigenvalue_Schnakenberg(n,m,epsilon):
    # Allocate memory for the eigenvalue that will be returned
    eigen_value = n*(n+1)
    # Calculate the linear term depending on the value of m
    if m == 0:
        next_term = ((4*(2*n+1))/(n*(n+1)))*(epsilon**2)
    elif:
        denom = (4**m)*np.fac(m)*np.fac(m-1)*np.fac(n-m)
        c = ((np.fac(m+n))/(denom))
        next_term+= -(2*n+1)*c**(epsilon**(2*m))
    # Add the next term in the asymptotic expansion
    eigen_value += next_term
    # Return the eigenvalue at hand
    return eigen_value
#------------------------------------------------------------------
# Function 4: "check_Turing_conditions_Scnakenberg"
# This function checks whether the provided parameters (a,b,gamma,d)
# satisfies the Turing conditions and in this case it returns the
# upper bound M and the lower bound L in the Turing analysis. 
#------------------------------------------------------------------
def check_Turing_conditions_Scnakenberg(a,b,d):
    # Calculate the four partial derivatives of
    # the Jacobian matrix used in the linear
    # stability analysis evaluated at the steady state
    f_u = ( (b - a) / (b + a) )
    f_v = ( (a + b)**2 )
    g_u = ( (-2*b) / (a+b) )
    g_v = -f_v
    # Calculate the determinant of the Jacobian matrix evaluated at the steady state
    det_A = (f_u*g_v) - (f_v*g_u)
    # Now calculate our four conditions
    cond_1 = (f_u+g_v<0)
    cond_2 = (det_A>0)
    cond_3 = ((d*f_u+g_v)>0)
    cond_4 = (((d*f_u+gv)**2-(4*d*det_A))>0)
    # Now we define the output depending on whether these conditions are met or not
    if cond_1 and cond_2 and cond_3 and cond_4:
        Turing_conditions = True
        discriminant = np.sqrt((d*f_u+g_v)**2-(4*d*det_A))
        L = ((d*f_u+g_v-discriminant)/(2*d))
        M = ((d*f_u+g_v+discriminant)/(2*d))        
    else:
        Turing_conditions = False
        L = 0
        M = 0
    # Finally, return these outputs
    return Turing_conditions,L,M
