# Implementação de um potencial de Lennard-Jones modificado para análise de interações entre partículas
O projeto tem como objetivo realizar uma simulação de Dinâmica Molecular, escrita em linguagem Python, implementando o potencial q-deformado, proposto por Alves, Farias e Santos (2022). 
O projeto está, no momento, em fase de documentação.

## Potencial q-deformado:
$$\Phi (\nu) = \frac{V_0 B_0}{(q-1)^2}\left[ \frac{1}{2} \nu^{-2(q-1)} - \nu^{-(q-1)} \right]\ (1)$$

É possível mostrar que este potencial (o qual será implementado na simulação) é equivalente ao 12-6 (para $q-1 = 2$). Da equação $(70)$ de [1], $\nu ^{-(1/3)}=(\frac{a}{r})$. Fazendo a substituição de $\nu$ e utilizando $q-1$ = 2:

$$\Phi (r) = \frac{V_0 B_0}{4}\left[ \frac{1}{2} \left(\frac{a}{r}\right)^{12} - \left(\frac{a}{r}\right)^{6} \right] = C \left[ \frac{1}{2} \left(\frac{a}{r}\right)^{12} - \left(\frac{a}{r}\right)^{6} \right] $$

Na distância de equilíbrio, $r = a$ e $\Phi (r) = -\epsilon$. Portanto,

$$\Phi(r) = C \left[ \frac{1}{2} \left(\frac{a}{a}\right)^{12} - \left( \frac{a}{a} \right)^6 \right] = C \left( \frac{1}{2} - 1\right) = -\frac{C}{2} = -\epsilon \implies C = 2 \epsilon$$
$$\therefore \Phi(r) = 2 \epsilon \left[ \frac{1}{2} \left(\frac{a}{r}\right)^{12} - \left(\frac{a}{r}\right)^{6} \right]$$

A uma distância $\sigma$, o potencial é nulo. Logo,

$$ 2 \epsilon \left[ \frac{1}{2} \left(\frac{a}{\sigma}\right)^{12} - \left(\frac{a}{\sigma}\right)^{6} \right] = 0 \implies \frac{1}{2} \left(\frac{a}{\sigma}\right)^{12} - \left(\frac{a}{\sigma}\right)^{6}  = 0 \implies \left( \frac{a}{\sigma} \right)^6 =  \frac{1}{2} \left(\frac{a}{\sigma}\right)^{12} \implies a^6 = 2 \sigma^6 \ , \ a^{12} = 4 \sigma ^{12}$$

$$ \therefore \Phi (r) = 2 \epsilon \left( \frac{1}{2} \frac{4 \sigma ^{12}}{r^{12}} - \frac{2 \sigma ^6}{r^6}\right) = 4 \epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right] $$

Mostrando que o potencial proposto é uma generalização do potencial de Lennard-Jones. Como foi identificado $C = \frac{V_0 B_0}{4} = 2 \epsilon$, e portanto, $4 \epsilon = 2 C$, e $q-1 = 2$, é possível escrever o potencial generalizado na seguinte forma:

$$\Phi (r) = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{2[{3(q-1)}]} - \left( \frac{\sigma}{r} \right)^{3(q-1)} \right] = 4 \epsilon \left[ \left( \frac{\sigma}{r} \right)^{2m} - \left( \frac{\sigma}{r} \right)^{m} \right] \ , \ m = 3(q-1)$$

Este será o potencial utilizado na simulação. A expressão para a força, portanto, será:

$$F= - \nabla \Phi (r) = 4 \epsilon \left(-2m \frac{\sigma ^{2m}}{r^{2m+1}}+ m \frac{\sigma ^m}{r^{m+1}}\right)$$

Referências bibliográficas: 
[1] ALVES, Á.S.; FARIAS, F.A.; SANTOS, R.N. 2022. A Proposal of an Isothermal q-EoS for Solids at High Pressures. Sitientibus Ser. Cienc. Fis. 18: 1-13.

[2] GLEIM, T. 2009. Velocity-Verlet-like algorithm for simulations of stochastic dynamics. arXiv preprint arXiv:0906.1863.

[3] JONES, J.E. 1924. On the determination of molecular fields. II. From the equation of state of a gas. Proc. R. Soc. A 106(738): 463-477.

[4] MADEIRA, L.; VITIELLO, S. 2012. A função radial de distribuição de pares para sistemas Lennard-Jones bidimensionais. Rev. Bras. Ensino Fís. 34(4).

[5] MARTÍNEZ, Leandro; BORIN, Ivana Aparecida; SKAF, Munir S. Fundamentos de Simulação por Dinâmica Molecular. In: MORGON, Nelson Henrique; COUTINHO, Kaline Rabelo (ed.). Métodos de química teórica e modelagem molecular. São Paulo: Editora Livraria da Física, 2007. cap. 12, p. 413-452.

[6] NAMBA, A. M.; SILVA, V. B. da; SILVA, C. H. T. P. da. Dinâmica molecular: teoria e aplicações em planejamento de fármacos. Eclética Química, Araraquara, v. 33, n. 4, p. 13-23, 2008. Disponível em: https://doi.org/10.26850/1678-4618eqj.v33.4.2008.p13-23. Acesso em: 14 out. 2026.

[7] PATERLINI, M.G.; FERGUSON, D.M. 1998. Constant temperature simulations using the Langevin equation with velocity Verlet integration. Chem. Phys. 236(1-3): 243-252.

[8] RAHMAN, A. 1964. Correlations in the Motion of Atoms in Liquid Argon. Phys. Rev. 136(2A): A405-A411.
[9] SCHWERDTFEGER, P.; WALES, D.J. 2024. 100 Years of the Lennard-Jones Potential. J. Chem. Theory Comput. 20(9): 3379-3405.
[10] SHEKAARI, A.; JAFARI, M. 2021. A concise introduction to molecular dynamics simulation: Theory and programming. arXiv preprint arXiv:2103.16944.
[11] TSALLIS, C. 1988. Possible generalization of Boltzmann-Gibbs statistics. J. Stat. Phys. 52: 479.

