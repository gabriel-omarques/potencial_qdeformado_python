import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MaxNLocator
from config import config

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['lines.linewidth'] = 1.8

def conversao(T_red, P_red, rho_red, cv_red, E_red):
    eps_kb = 119.8
    sigma = 3.405
    M = 39.95
    R = 8.31446
    T_real = T_red * eps_kb
    rho_real = rho_red * (M / (0.602214 * sigma**3))
    P_real = P_red * ((eps_kb * 138.0649) / (sigma**3))
    cv_real = (cv_red + 1.5) * R
    U_real = E_red * eps_kb * R / 1000.0
    return T_real, P_real, rho_real, cv_real, U_real

def block_average(dados, n_blocos=20):
    dados = np.asarray(dados)
    n = len(dados)
    if n < 2:
        return (dados[0] if n == 1 else 0.0), 0.0
    n_blocos = max(1, min(n_blocos, n))
    blocos = np.array_split(dados, n_blocos)
    medias_bloco = np.array([b.mean() for b in blocos if len(b) > 0])
    media_global = dados.mean()
    erro = medias_bloco.std(ddof=1) / np.sqrt(len(medias_bloco)) if len(medias_bloco) > 1 else 0.0
    return media_global, erro

def block_cv(pot_series, N, T0, n_blocos=20):
    pot_series = np.asarray(pot_series)
    n = len(pot_series)
    if n < 4:
        return 0.0, 0.0
    n_blocos = max(2, min(n_blocos, n // 2))
    blocos = np.array_split(pot_series, n_blocos)
    cvs = np.array([
        (b.var() * (N**2)) / (N * (T0**2))
        for b in blocos if len(b) > 1
    ])
    if len(cvs) > 1:
        return cvs.mean(), cvs.std(ddof=1) / np.sqrt(len(cvs))
    return (cvs[0] if len(cvs) == 1 else 0.0), 0.0

def _rdf_histograma_frame(pos, aresta, bins, N):
    delta = pos[:, None, :] - pos[None, :, :]
    delta -= aresta * np.round(delta / aresta)
    r = np.sqrt(np.sum(delta * delta, axis=-1))
    iu = np.triu_indices(N, k=1)
    r_pairs = r[iu]
    rmax = bins[-1]
    r_pairs = r_pairs[r_pairs < rmax]
    hist, _ = np.histogram(r_pairs, bins=bins)
    return hist * 2

def calc_analise():
    resultado = {
        'C_v': 0.0, 'C_v_err': 0.0, 'T_med': 0.0, 'T_err': 0.0,
        'P_med': 0.0, 'P_err': 0.0, 'K_med': 0.0, 'K_err': 0.0,
        'U_med': 0.0, 'U_err': 0.0, 'E_med': 0.0, 'E_err': 0.0,
        'dados_reais': conversao(0.0, 0.0, config.rho, 0.0, 0.0),
        'r': np.array([0]), 'gr': np.array([0]),
        'tpz': None, 'kpe': None,
        'u_tail': config.u_tail, 'p_tail': config.p_tail,
    }

    try:
        tpz = np.loadtxt(os.path.join(BASE_DIR, 'tpz.out'))
        kpe = np.loadtxt(os.path.join(BASE_DIR, 'kpe.out'))
        resultado['tpz'], resultado['kpe'] = tpz, kpe

        ign = config.ign
        T_series = tpz[ign:, 1]
        P_series_bruto = tpz[ign:, 2]
        K_series = kpe[ign:, 1]
        U_series_bruto = kpe[ign:, 2]
        E_series_bruto = kpe[ign:, 3]

        T_med, T_err = block_average(T_series)
        P_med_bruto, P_err = block_average(P_series_bruto)
        K_med, K_err = block_average(K_series)
        U_med_bruto, U_err = block_average(U_series_bruto)
        E_med_bruto, E_err = block_average(E_series_bruto)

        P_med = P_med_bruto + config.p_tail
        U_med = U_med_bruto + config.u_tail
        E_med = E_med_bruto + config.u_tail

        C_v, C_v_err = block_cv(U_series_bruto, config.N, config.T_0)

        resultado.update({
            'C_v': C_v, 'C_v_err': C_v_err, 'T_med': T_med, 'T_err': T_err,
            'P_med': P_med, 'P_err': P_err, 'K_med': K_med, 'K_err': K_err,
            'U_med': U_med, 'U_err': U_err, 'E_med': E_med, 'E_err': E_err,
        })
        resultado['dados_reais'] = conversao(T_med, P_med, config.rho, C_v, E_med)

    except (OSError, IndexError, ValueError) as e:
        print(f"[calc_analise] Aviso: Erro ao ler logs ({e}).")
        return resultado

    try:
        with open(os.path.join(BASE_DIR, 'pos.xyz'), 'r') as f:
            linhas = f.readlines()
    except OSError as e:
        print(f"[calc_analise] Aviso: nao foi possivel ler pos.xyz ({e}).")
        return resultado

    passos = []
    idx = 0
    while idx < len(linhas):
        n_atoms = int(linhas[idx])
        idx += 2
        coords = []
        for _ in range(n_atoms):
            partes = linhas[idx].split()
            c = partes[1:] if partes[0].isalpha() else partes
            coords.append([float(c[0]), float(c[1]), float(c[2])])
            idx += 1
        passos.append(np.array(coords))

    fotos_ignoradas = config.ign // 100
    trajetoria = np.array(passos)[fotos_ignoradas:]
    n_frames = len(trajetoria)

    if n_frames <= 1:
        return resultado

    dr = 0.1
    rmax = config.aresta / 2.0
    n_bins = int(round(rmax / dr))
    bins = np.linspace(0, n_bins * dr, n_bins + 1)
    hist_total = np.zeros(n_bins)
    
    for t in range(n_frames):
        hist_total += _rdf_histograma_frame(trajetoria[t], config.aresta, bins, config.N)

    r_centros = bins[:-1] + dr / 2.0
    g_r = np.zeros_like(hist_total)
    for i, r in enumerate(r_centros):
        vol_casca = 4.0 * np.pi * r**2 * dr
        g_r[i] = (hist_total[i] / n_frames) / (config.N * config.rho * vol_casca)

    resultado['r'] = r_centros[:-1]
    resultado['gr'] = g_r[:-1]
    return resultado

def graficos():
    resultado = calc_analise()
    T_real, P_real, rho_real, cv_real, U_real = resultado['dados_reais']

    P_nist_mpa = 10.023 
    U_nist = 3.1270      
    Cv_nist = 13.146      

    P_real_mpa = P_real / 10.0

    erro_P = abs(P_real_mpa - P_nist_mpa) / P_nist_mpa * 100.0
    erro_U = abs(U_real - U_nist) / abs(U_nist) * 100.0  
    erro_Cv = abs(cv_real - Cv_nist) / Cv_nist * 100.0

    print(f"\n{'='*65}")
    print(f"      Comparacao com os dados do NIST (Gas)")
    print(f"{'='*65}")
    print(f" Propriedade         | Simulacao      | NIST           | Erro (%)")
    print(f"-----------------------------------------------------------------")
    print(f" Pressao (MPa)       | {P_real_mpa:14.4f} | {P_nist_mpa:14.4f} | {erro_P:7.2f}%")
    print(f" Energia (kJ/mol)    | {U_real:14.4f} | {U_nist:14.4f} | {erro_U:7.2f}%")
    print(f" C_v (J/mol.K)       | {cv_real:14.4f} | {Cv_nist:14.4f} | {erro_Cv:7.2f}%")
    print(f"{'='*65}")

    print(f"{'-'*65}")
    print(f" Medias em unidades reduzidas")
    print(f"   T*  = {resultado['T_med']:.4f} ")
    print(f"   P*  = {resultado['P_med']:.4f} ")
    print(f"   U*  = {resultado['U_med']:.4f} ")
    print(f"   E*  = {resultado['E_med']:.4f} ")
    print(f"   C_v = {resultado['C_v']:.4f} ")
    print(f"{'='*65}\n")

    tpz, kpe = resultado['tpz'], resultado['kpe']
    r, gr = resultado['r'], resultado['gr']

    if tpz is not None and kpe is not None:
        fig2, ax2 = plt.subplots(figsize=(6, 4.5))
        
        r_teo = np.linspace(0.85, max(r), 500)
        sr_teo = config.sigma / r_teo
        srm_teo = sr_teo**config.m
        u_teo = 4.0 * config.epsilon * (srm_teo**2 - srm_teo)
        
        u_teo = np.clip(u_teo, a_min=-np.inf, a_max=700.0)
        g_teo = np.exp(-u_teo / config.T_0)
        

        ax2.plot(r_teo, g_teo, color='black', linestyle='-', linewidth=1.5, alpha=0.7, label='$ g(r) \approx \\exp(-u/T_0) $', zorder=1)
        
        ax2.plot(r, gr, marker='+', linestyle='None', color='red', markersize=5, label='Simulação', zorder=2)
        
        ax2.axhline(1.0, color='gray', ls='--', lw=1.2)
        ax2.set_xlabel(r'$r/\sigma$')
        ax2.set_ylabel(r'$g(r)$')
        ax2.set_xlim(0, max(r))
        ax2.set_ylim(0, max(gr) * 1.1)
        
        ax2.legend(loc='upper right', frameon=False)
        
        ax2.xaxis.set_major_locator(MaxNLocator(prune='lower'))
        ax2.xaxis.set_minor_locator(AutoMinorLocator())
        ax2.yaxis.set_minor_locator(AutoMinorLocator())
        plt.tight_layout()
        plt.savefig(os.path.join(BASE_DIR, 'g_r.png'), dpi=300)

        plt.show()