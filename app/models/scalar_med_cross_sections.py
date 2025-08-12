# ↑ Copy from the last output

import numpy as np


@np.vectorize
def gamma_scalar(g_s_psi, m_f, M_med):
    return g_s_psi**2*np.sqrt(1 - 4*m_f**2/M_med**2)/(np.pi*M_med)

@np.vectorize
def gamma_fermion(g_s_psi, m_f, M_med):
    return (1/16)*g_s_psi**2*np.sqrt(1 - 4*m_f**2/M_med**2)*(M_med**2 - 4*m_f**2)/(np.pi*M_med)

@np.vectorize
def gamma_vector(g_s_psi, m_f, M_med):
    return (1/64)*g_s_psi**2*np.sqrt(1 - 4*m_f**2/M_med**2)*(M_med**4 - 8*M_med**2*m_f**2 + 16*m_f**4)/(np.pi*M_med*m_f**4)


gamma_func_dict = {'scalar': gamma_scalar,
                   'fermion': gamma_fermion,
                   'vector': gamma_vector,
}
## 

@np.vectorize
def scalar_DM_scalar_med(M_med, g_ps_psi, g_ps_chi, s, g_s_chi, g_s_psi, m_f, Gamma_func):
    
    g_ps_psi = g_s_psi # Equal couplings for scalar and pseudo-scalar mediators
    g_ps_chi = g_s_chi # Equal couplings for scalar and pseudo-scalar mediators
    Bi = Gamma_func(g_ps_psi, m_f, M_med) # Initial states branching ratio
    Bf = Gamma_func(g_ps_chi, m_f, M_med) # Final states branching ratio
    Gamma = Bi + Bf

    return 16*g_s_chi**2*(g_ps_psi**2 + g_s_psi**2)*np.sqrt(-4*m_f**2 + s)/(np.sqrt(s)*(Gamma**2*M_med**2 + M_med**4 - 2*M_med**2*s + s**2))

@np.vectorize
def fermion_DM_scalar_med(M_med, g_ps_psi, g_ps_chi, s, g_s_chi, g_s_psi, m_f, Gamma_func):
    
    g_ps_psi = g_s_psi # Equal couplings for scalar and pseudo-scalar mediators
    g_ps_chi = g_s_chi # Equal couplings for scalar and pseudo-scalar mediators
    Bi = Gamma_func(g_ps_psi, m_f, M_med) # Initial states branching ratio
    Bf = Gamma_func(g_ps_chi, m_f, M_med) # Final states branching ratio
    Gamma = Bi + Bf

    
    return (1/2)*np.sqrt(-4*m_f**2 + s)*(-4*g_s_chi**2*m_f**2*(g_ps_psi**2 + g_s_psi**2) + s*(g_ps_chi**2*g_ps_psi**2 + g_ps_chi**2*g_s_psi**2 + g_ps_psi**2*g_s_chi**2 + g_s_chi**2*g_s_psi**2))/(np.sqrt(s)*(Gamma**2*M_med**2 + M_med**4 - 2*M_med**2*s + s**2))

@np.vectorize
def vector_DM_scalar_med(M_med, g_ps_psi, g_ps_chi, s, g_s_chi, g_s_psi, m_f, Gamma_func):

    g_ps_psi = g_s_psi # Equal couplings for scalar and pseudo-scalar mediators
    g_ps_chi = g_s_chi # Equal couplings for scalar and pseudo-scalar mediators
    Bi = Gamma_func(g_ps_psi, m_f, M_med) # Initial states branching ratio
    Bf = Gamma_func(g_ps_chi, m_f, M_med) # Final states branching ratio
    Gamma = Bi + Bf

    
    return (3/4)*g_s_chi**2*(g_ps_psi**2 + g_s_psi**2)*np.sqrt(-4*m_f**2 + s)/(np.sqrt(s)*(Gamma**2*M_med**2 + M_med**4 - 2*M_med**2*s + s**2))


xsec_func_dict = {'scalar': scalar_DM_scalar_med,
                   'fermion': fermion_DM_scalar_med,
                   'vector': vector_DM_scalar_med,
}