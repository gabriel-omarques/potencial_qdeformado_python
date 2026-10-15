import numpy as np

class Config:
    SEED = 42 
    DIM = 3 
    N = 1000
    rho = 0.1        
    vol = N / rho 
    aresta = vol**(1/3) 
    NS = 1000000
    dt = 0.001 
    ign = 10000       
    T_0 = 2.5        
    tau = 0.01 
    epsilon = 1.0
    sigma = 1.0
    r_ct = 2.5 * sigma 
    sr_ct = sigma / r_ct
    q = 3.0 
    m = 3.0 * (q - 1.0)
    sr_ct_m = sr_ct**m
    sr_ct_2m = sr_ct_m**2
    u_ct = 4.0 * epsilon * (sr_ct_2m - sr_ct_m)
    
    u_tail = 8.0 * np.pi * rho * epsilon * (sigma**3) * ((1.0 / (2.0 * m - 3.0)) * sr_ct**(2.0 * m - 3.0) - (1.0 / (m - 3.0)) * sr_ct**(m - 3.0))
    p_tail = (8.0 / 3.0) * np.pi * (rho**2) * epsilon * m * (sigma**3) * ((2.0 / (2.0 * m - 3.0)) * sr_ct**(2.0 * m - 3.0) - (1.0 / (m - 3.0)) * sr_ct**(m - 3.0))

config = Config()