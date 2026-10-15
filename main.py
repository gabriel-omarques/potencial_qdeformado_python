import time
import os
import numpy as np
from config import config
from motor import forcas, estatisticas
from analise import graficos

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def main():
    np.random.seed(config.SEED)
    pos = np.zeros((config.DIM, config.N))
    n_lado = int(np.ceil(config.N**(1/3)))
    espacamento = 1.0 / n_lado

    idx = 0
    for ix in range(n_lado):
        for iy in range(n_lado):
            for iz in range(n_lado):
                if idx < config.N:
                    pos[0, idx] = -0.5 + (ix + 0.5) * espacamento
                    pos[1, idx] = -0.5 + (iy + 0.5) * espacamento
                    pos[2, idx] = -0.5 + (iz + 0.5) * espacamento
                    idx += 1

    vel = np.random.uniform(-0.5, 0.5, (config.DIM, config.N))
    vel -= np.mean(vel, axis=1, keepdims=True)

    print("Inicio da simulacao.")
    start_time = time.time()

    tpz_buf = np.zeros((config.NS, 3))
    kpe_buf = np.zeros((config.NS, 4))

    c1 = np.exp(-config.dt / config.tau)
    c2 = np.sqrt((1.0 - c1**2) * config.T_0) / config.aresta 

    caminho_pos = os.path.join(BASE_DIR, 'pos.xyz')

    with open(caminho_pos, 'w') as f_xyz:
        pos[pos > 0.5] -= 1.0
        pos[pos < -0.5] += 1.0

        acc, pot, vrl = forcas(pos, config.N, config.DIM, config.aresta,
                               config.r_ct, config.epsilon, config.sigma, config.u_ct, config.m)

        for step in range(1, config.NS + 1):
            pos += config.dt * vel + 0.5 * acc * config.dt**2

            ruido = np.random.normal(0.0, 1.0, (config.DIM, config.N))
            ruido -= np.mean(ruido, axis=1, keepdims=True)
            
            vel_termo = c1 * vel + c2 * ruido
            vel_meia = vel_termo + 0.5 * config.dt * acc

            pos[pos > 0.5] -= 1.0
            pos[pos < -0.5] += 1.0

            acc, pot, vrl = forcas(pos, config.N, config.DIM, config.aresta,
                                   config.r_ct, config.epsilon, config.sigma, config.u_ct, config.m)
            
            vel = vel_meia + 0.5 * config.dt * acc

            k_avg, T_i = estatisticas(vel, config.aresta, config.DIM, config.N)
            
            p_avg = np.mean(pot)
            pres = config.rho * T_i + vrl / config.vol

            tpz_buf[step - 1] = (step * config.dt, T_i, pres)
            kpe_buf[step - 1] = (step * config.dt, k_avg, p_avg, k_avg + p_avg)

            if step % 2000 == 0:
                f_xyz.write(f"{config.N}\nPasso {step}\n")
                np.savetxt(f_xyz, (pos * config.aresta).T, fmt='Ar %.4f %.4f %.4f')

            if step % (config.NS // 100) == 0:
                print(f"{(step / config.NS) * 100:.1f}%")

    np.savetxt(os.path.join(BASE_DIR, 'tpz.out'), tpz_buf, fmt='%.6e')
    np.savetxt(os.path.join(BASE_DIR, 'kpe.out'), kpe_buf, fmt='%.6e')

    print(f"Tempo de simulacao: {time.time() - start_time:.2f} s")
    graficos()

if __name__ == "__main__":
    main()