import numpy as np
from math import copysign
from numba import njit, prange

@njit(fastmath=True, parallel=True)
def forcas(pos, N, dim, aresta, r_ct, epsilon, sigma, u_ct, m):
    acc = np.zeros((dim, N))
    pot = np.zeros(N)
    vrl_por_particula = np.zeros(N)
    
    sigma2 = sigma * sigma
    r_ct2 = r_ct * r_ct
    m_metade = m / 2.0

    for i in prange(N): 
        for j in range(N): 
            if i != j: 
                r2_ij = 0.0 
                R_ij = np.zeros(dim) 

                for k in range(dim):
                    R_ij[k] = pos[k, i] - pos[k, j] 
                    if abs(R_ij[k]) > 0.5:
                        R_ij[k] -= copysign(1.0, R_ij[k]) 
                    r_k = aresta * R_ij[k] 
                    r2_ij += r_k * r_k 

                if r2_ij < r_ct2: 
                    sr2 = sigma2 / r2_ij 
                    srm = sr2 ** m_metade
                    sr2m = srm * srm

                    u = 4.0 * epsilon * (sr2m - srm) - u_ct 
                    pot[i] += 0.5 * u 

                    du = (4.0 * epsilon * m / r2_ij) * (2.0 * sr2m - srm) 

                    vrl_por_particula[i] -= 0.5 * du * r2_ij 

                    for k in range(dim):
                        acc[k, i] += du * R_ij[k] 

    vrl = vrl_por_particula.sum() 
    return acc, pot, -vrl / dim

def estatisticas(vel, aresta, dim, N):
    v2_total = np.sum((vel * aresta)**2)
    k_avg = 0.5 * v2_total / N
    N_f = dim * N - dim         
    temp_i = v2_total / N_f
    return k_avg, temp_i