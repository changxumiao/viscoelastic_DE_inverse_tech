# Introduction
A priori knowledge on micro-level scale kinetics facilitate the search for Prony Series (c.f. [discussion][search_PS_from_kinetics]). However, this knowledge is usually unavailable. Searching for Prony Series from measurement is an alternative to solve the problem. 

We categorize the measurement into two types:
1. Standard test, such as uniaxial tension, pure shear, etc.
2. Non-standard test, such as force-strain data acquired in working condition.

# Uniaxial tension
From article [^physical_gel], the nominal stress is
[^physical_gel]: Hui, C.-Y., Cui, F., Zehnder, A., & Vernerey, F. J. (2021). Physically motivated models of polymer networks with dynamic cross-links: Comparative study and future outlook. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 477(2255), 20210608. https://doi.org/10.1098/rspa.2021.0608

$$ P=\mu \{[\rho_{p}+n(t)][\lambda(t)-\frac{1}{\lambda^2(t)}]+\gamma_\infty\int_0^t\phi_B(t-\tau)[\frac{\lambda(t)}{\lambda^2(\tau)}-\frac{\lambda(\tau)}{\lambda^2(t)}]d\tau\}$$

With time series of nominal stress (P) and stretch ratio ($\lambda$) available.

Assume
$$\phi_B(t-\tau)\approx\sum_n g^{(n)}e^{-\frac{t-\tau}{\tau^{(n)}}}$$

$$n(t)=\int_{-\infty}^0\phi_B(t-\tau)d\tau\approx \sum_n g^{(n)}\tau^{(n)}e^{-\frac{t}{\tau^{(n)}}}$$

The question becomes deconvolution: finding lest number of $g^{(n)}$ that satisfy the above equation.

# Generate signal
| Parameter | Value |
| --------- | ----- |
| $\mu\rho_p(\rm{kPa})$ | $2.347$ |
| $\mu\gamma_\infty (\rm{kPas^{-1}})$ | $19.91$ |
| $\alpha$ | $1.634$ |
| $t_B(\rm{s})$ | $0.461$ |

```matlab
rhop0mu = 2.347e3;  %[Pa]
gammarhod0mu = 19.91e3; %[Pa]
tB = 0.4618;    %[s]
alpha =  1.634;
```
## Row data from Power Law
```matlab
time_load1 = 0.3/0.003;%    [s]
t_max = 2*time_load1;%   [s]
time_series = linspace(0,t_max,1000);
dt = t_max/length(time_series);
lambda_t = 1+0.003*time_series.*(time_series<time_load1)+((-0.003)*(time_series-time_load1)+0.003*time_load1).*(time_series>=time_load1);%    [1]

an1 = (1-(1-alpha)*time_series/tB).^(1/(1-alpha));
an1_int = tB/(2-alpha)*(1+(alpha-1)*time_series/tB).^((2-alpha)/(1-alpha));

term1 = rhop0mu*(lambda_t-1./lambda_t.^2);
term2 = gammarhod0mu*an1_int.*(lambda_t-1./lambda_t.^2);
temp_cache = conv(an1,1./lambda_t.^2,'full')*dt;
term3 = gammarhod0mu*lambda_t.*cache(1:length(time_series));
cache = conv(an1,lambda_t,'full')*dt;
term4 = -gammarhod0mu*1./lambda_t.^2.*cache(1:length(time_series));
P = term1+term2+term3+term4;% Nominal stress
s = P.*lambda_t.^2% Cauchy stress
```
The above can be replaced by the alternatives of power law: Prony Series
```matlab
g1 = 0.9928;
tau1 = 0.6424;
an2 = g1*exp(-time_series/tau1);    %Prony Series alternative
```
## Interpolation
```matlab
P_rep = interp1(lambda_t,P,xq)
```
# Matrix form of convolution
The general form of convolution is

$$f=h*u$$
The $f$ is output signal, $u$ is input signal, $h$ is the filter (impulse response).

```matlab
function h_matrix1 = conv_matrix(h,N)
% return a square matrix of convolution operator, with size equal to N
M = length(h);
L = M+N-1;
h_aug = [h,zeros(1,L-M)];
h_matrix = zeros(L,L);
for ii=1:L
    h_matrix(:,ii) = circshift(h_aug,ii-1)';
end
h_matrix1 = h_matrix(1:N,1:N);
end
```
Example
```matlab
u = [1 0 2];
h = [2 7];
f = conv_matrix(h,length(u))*u'
```

# Deconvolution from stress-strain data

$$ P= \mu\rho_{p}[\lambda(t)-\frac{1}{\lambda^2(t)}]+\mu\gamma_\infty[\lambda(t)-\frac{1}{\lambda^2(t)}]\sum_n g^{(n)}\tau^{(n)}e^{-\frac{t}{\tau^{(n)}}}+\lambda(t) \mu\gamma_\infty\int_0^t\sum_n g^{(n)}e^{-\frac{t-\tau}{\tau^{(n)}}}\frac{1}{\lambda^2(\tau)}d\tau-\frac{1}{\lambda^2(t)}\mu\gamma_\infty\int_0^t\sum_n g^{(n)}e^{-\frac{t-\tau}{\tau^{(n)}}}\lambda(\tau)d\tau$$

Its abstract form is

$$P=f*u+b$$

Impose a dictionary, the form is

$$P=Hg+b$$

to find $g$ with accuracy and sparsity.
```matlab
b = rhop0mu*(lambda_t-1./lambda_t.^2);
H = gammarhod0mu*(lambda_t-1./lambda_t.^2)'.*tau.*A+gammarhod0mu*lambda_t'.*conv_matrix(1./lambda_t.^2,length(time_series))*dt*A-gammarhod0mu*(1./lambda_t.^2)'.*conv_matrix(lambda_t,length(time_series))*dt*A;
signal = P'-b';
```
Create dictionary, normalize-le and apply matching pursuit
```matlab
%% Normalized H
for ii=1:size(H,2)
    norm_H2(ii) = norm(H(:,ii),2);
    H_norm(:,ii) = H(:,ii)/norm_H2(ii);
end
g = matchingpursuit(signal,H_norm);
G_recover = g./norm_H2';
tau_recover = tau(find(G_recover));
```
# Result



[search_PS_from_kinetics]: ./search_PS_from_kinetics.md