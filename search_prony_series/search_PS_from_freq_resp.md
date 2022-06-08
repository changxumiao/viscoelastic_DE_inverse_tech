# Introduction
Characterizing viscosity via frequency response is another way to find Prony series.The doc is committed to search Prony Series from measured frequency response.[^Control][^cell]

We denote:
- shear state of material by $\gamma$
- time-dependent shear modulus by $G(t)$
- shear stress by $s(t)$

Within the scope of linear viscoelasitc material, we have

$$s(t)=\int_0^tG(t-\tau)\frac{ds}{d\tau} d\tau$$

# General Maxwell Model
The generalized Maxwell model consists of K parallel spring-damper elements (viscous part) and a supplementary parallel linear elastic spring (elastic part).[^Kruas]

$$s(t)=\int_0^t[G_\infty +\sum_{k=1}^KG_k \rm{exp}(-\frac{t-\tau}{\tau_k})]\frac{ds}{d\tau} d\tau$$

where

$$\tau_k=\frac{\eta_k}{G_k}$$

The function 

$$G(t-\tau)=G_\infty +\sum_{k=1}^K G_k \rm{exp}(-\frac{t-\tau}{\tau_k})$$

is called Prony Series

If we assume network kinetics is **independent to deformation**, we can have
$$G_\infty=G_0$$

# Frequency domain

The time dependent relaxation properties are viewed in frequency domain: complex moduli.

The Fourier transformation[^wiki_FT] of Prony series goes as

$$G^*(\omega)=G_\infty+\sum_{k=1}^K G_k \frac{\omega^2\tau^2_k}{1+\omega^2\tau^2_k}+i\sum_{k=1}^K G_k \frac{\omega^2\tau_k}{1+\omega^2\tau^2_k}$$

Namely

$$G^*(\omega)=G_0(1+\sum_{k=1}^K g_k \frac{\omega^2\tau^2_k}{1+\omega^2\tau^2_k}+i\sum_{k=1}^K g_k \frac{\omega^2\tau_k}{1+\omega^2\tau^2_k})$$

$G'(\omega)$, the real part of $G^*(\omega)$, is called the *storage modulus*; $G''(\omega)$, the imaginary part of $G^*(\omega)$, is called the *loss modulus*.

# Abstract problem
Find least number of $\{g_k,\tau_k\}$, subject to

$$G' = G_0(1+\sum_{k=1}^K g_k \frac{\omega^2\tau^2_k}{1+\omega^2\tau^2_k})$$

# Generating dictionary
```matlab
omega_f = logspace(-6,4,1000);
tau = logspace(-2,4,1000);
A = zeros(length(omega_f),length(tau));
for ii = 1:length(tau)
    A(:,ii) = (omega_f.^2.*tau(ii)^2)./(1+omega_f.^2.*tau(ii)^2)';
    norm_2(ii) = norm(A(:,ii),2);
    A(:,ii) = A(:,ii)/norm_2(ii);
end
```

# Searchg PS by Matching pursuit
```matlab
G_real = 
signal = G_real-G_0;
g = matchingpursuit(signal,A);
G_recover = g./norm_2';
tau_recover = tau(find(G_recover));
```


[^Kruas]: Kraus, M. A., Schuster, M., Kuntsche, J., Siebert, G., & Schneider, J. (2017). Parameter identification methods for visco- and hyperelastic material models. Glass Structures & Engineering, 2(2), 147–167. https://doi.org/10.1007/s40940-017-0042-9

[^wiki_FT]: https://en.wikipedia.org/wiki/Fourier_transform#:~:text=%3D%200.-,Distributions%2C%20one%2Ddimensional,-%5Bedit%5D

[^Control]: Grindy, S. C., Learsch, R., Mozhdehi, D., Cheng, J., Barrett, D. G., Guan, Z., Messersmith, P. B., & Holten-Andersen, N. (2015). Control of hierarchical polymer mechanics with bioinspired metal-coordination dynamics. Nature Materials, 14(12), 1210–1216. https://doi.org/10.1038/nmat4401

[^cell]: Hang, J.-T., Xu, G.-K., & Gao, H. (2022). Frequency-dependent transition in power-law rheological behavior of living cells. Science Advances, 8(18), eabn6093. https://doi.org/10.1126/sciadv.abn6093
