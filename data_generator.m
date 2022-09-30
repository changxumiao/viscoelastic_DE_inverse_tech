clear all

%loading type
unaxial_tension = 0;
pure_shear = 1;

markerinterval = 200;
check = 1;
back2real = 1;


mu = 12.42e3;
N1 = 0.1234;%fraction of permenant bonds
N2= (1-N1);

%% Prony Series
tau1 = 0.1;
tau2 = 4.2;
tau3 = 33.78;
tau4 = 3134;

g1 = 0.1234;
g2 = 0.1234;
g3 = 0.1234;
g4 = 0.1234;

tau1_dimless = tau1/tau1;
tau2_dimless = tau2/tau1;
tau3_dimless = tau3/tau1;
tau4_dimless = tau4/tau1;

%% loading info
disp_max_dimless = 1;
loading_rate = 0.03;
loading_rate_dimless = loading_rate*tau1;
loading_time = disp_max_dimless/loading_rate_dimless;
deloading_time = loading_time;
time_sum = loading_time+deloading_time;
%%
time_series = linspace(0,time_sum,9999);
dt = time_sum/length(time_series);

FxX = ones(1,length(time_series));
FxY = zeros(1,length(time_series));
FxZ = zeros(1,length(time_series));
FyX = zeros(1,length(time_series));
FyY = ones(1,length(time_series));
FyZ = zeros(1,length(time_series));
FzX = zeros(1,length(time_series));
FzY = zeros(1,length(time_series));
FzZ = ones(1,length(time_series));
F = zeros(3,3,length(time_series));
H = zeros(3,3,length(time_series));
%% loading type
if unaxial_tension
    FxX = 1 + loading_rate*time_series.*(time_series<loading_time)+((-loading_rate)*(time_series-loading_time)+loading_rate*loading_time).*(time_series>=loading_time);%    [1]
    FyY = 1./FxX.^0.5;
    FzZ = 1./FxX.^0.5;
end

if pure_shear
    FxY = (loading_rate*time_series.*(time_series<loading_time)+((-loading_rate)*(time_series-loading_time)+loading_rate*loading_time).*(time_series>=loading_time));%    [1]
end

%% Assemble matrix of deformation gradient
for ii = 1:1:length(time_series)
    F(:,:,ii) = [FxX(ii) FxY(ii) FxZ(ii); FyX(ii) FyY(ii) FyZ(ii); FzX(ii) FzY(ii) FzZ(ii)];
    H(:,:,ii) = inv(F(:,:,ii));
    H_HT(:,:,ii) = H(:,:,ii)*H(:,:,ii)';
    J(ii) = det(F(:,:,ii));
end

if check
    figure
    plot(time_series,squeeze(F(1,1,:)),'o','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(1,2,:)),'+','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(1,3,:)),'*','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(2,1,:)),'.','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(2,2,:)),'x','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(2,3,:)),'s','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(3,1,:)),'d','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(3,2,:)),'^','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(F(3,3,:)),'v','MarkerIndices',1:markerinterval:length(time_series));hold on;
    legend('FxX','FxY','FxZ','FyX','FyY','FyZ','FzX','FzY','FzZ')
    xlabel('time')
    ylabel('Deformation gradient F')
    title('Deformation gradient F')
end
if check
    figure
    plot(time_series,squeeze(H(1,1,:)),'o','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(1,2,:)),'+','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(1,3,:)),'*','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(2,1,:)),'.','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(2,2,:)),'x','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(2,3,:)),'s','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(3,1,:)),'d','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(3,2,:)),'^','MarkerIndices',1:markerinterval:length(time_series));hold on;
    plot(time_series,squeeze(H(3,3,:)),'v','MarkerIndices',1:markerinterval:length(time_series));hold on;
    legend('HxX','HxY','HxZ','HyX','HyY','HyZ','HzX','HzY','HzZ')
    xlabel('time')
    ylabel('inverse of deformation gradient (H)')
    title('inverse of deformation gradient (H)')
end
%%
an1 = g1/tau1_dimless*exp(-time_series/tau1_dimless)+g2/tau2_dimless*exp(-time_series/tau2_dimless)+g3/tau3_dimless*exp(-time_series/tau3_dimless)+g4/tau4_dimless*exp(-time_series/tau4_dimless);
an1_int = g1*exp(-time_series/tau1_dimless)+g2*exp(-time_series/tau2_dimless)+g3*exp(-time_series/tau3_dimless)+g4*exp(-time_series/tau4_dimless);
if check
    figure
    plot(time_series,an1)
    xlabel('time')
    ylabel('an1')
    title('Prony Series')
end

if check
    figure
    plot(time_series,an1_int)
    xlabel('time')
    ylabel('an1 int')
    title('Integration of Prony Series on \tau')
end
    
%% Elastic contribution
for ii = 1:1:length(time_series)
    S_elastic(:,:,ii) = N1*eye(3,3);
end
%% Viscous contribution (1)

for ii = 1:1:length(time_series)
    S_viscous1(:,:,ii) = N2*eye(3,3)*an1_int(ii);
end

%% Viscous contribution (2)
for jj = 1:1:3
    for kk = jj:1:3
        cache = conv(an1,squeeze(H_HT(jj,kk,:))','full')*dt;
        S_viscous2(jj,kk,:) = N2*cache(1:length(time_series));
    end
end

%% sum of elastic and viscous contributions
S = S_elastic + S_viscous1 + S_viscous2;

for ii = 1:1:length(time_series)
    P(:,:,ii) = F(:,:,ii)*S(:,:,ii);
    s(:,:,ii) = F(:,:,ii)*S(:,:,ii)*F(:,:,ii)'/J(ii);
end

%% deviatoric stress
for ii = 1:1:length(time_series)
    p_auxiliary(ii) = (1/3)*trace(s(:,:,ii));
    s_dev(:,:,ii) = s(:,:,ii) - p_auxiliary(ii)*eye(3,3);
    P_dev(:,:,ii) = J(ii)*s_dev(:,:,ii)*H(:,:,ii)';
    S_dev(:,:,ii) = J(ii)*H(:,:,ii)*s_dev(:,:,ii)*H(:,:,ii)';
end


%%
if check
    figure
    plot(time_series,J)
    xlabel('time')
    ylabel('J')
    title('Volumetric change')

    figure
    plot(time_series,p_auxiliary)
    xlabel('time')
    ylabel('p_{auxil}')
    title('hydrostatic pressure')
end


%%
figure
plot(time_series,squeeze(P_dev(1,1,:)),'o','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(1,2,:)),'+','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(1,3,:)),'*','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(2,1,:)),'.','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(2,2,:)),'x','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(2,3,:)),'s','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(3,1,:)),'d','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(3,2,:)),'^','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(P_dev(3,3,:)),'v','MarkerIndices',1:markerinterval:length(time_series));hold on;
xlabel('time')
ylabel('Nominal stress P (deviatoric)')
legend('P_d11','P_d12','P_d13','P_d21','P_d22','P_d23','P_d31','P_d32','P_d33')

figure
plot(time_series,squeeze(s_dev(1,1,:)),'o','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(s_dev(1,2,:)),'+','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(s_dev(1,3,:)),'*','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(s_dev(2,2,:)),'.','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(s_dev(2,3,:)),'x','MarkerIndices',1:markerinterval:length(time_series));hold on;
plot(time_series,squeeze(s_dev(3,3,:)),'s','MarkerIndices',1:markerinterval:length(time_series));hold on;
xlabel('time')
ylabel('Cauchy Stress s')
title('deviatoric Cauchy stress')
legend('s_d11','s_d12','s_d13','s_d22','s_d23','s_d33')

if unaxial_tension
    figure
    plot(FxX,squeeze(P_dev(1,1,:)),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('FxX')
    ylabel('P_d11')
end

if pure_shear
    figure
    plot(FxY,squeeze(s_dev(1,2,:)),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('FxY')
    ylabel('s_d12')
end

%% back to reality

if unaxial_tension
    figure
    plot(time_series*tau1,squeeze(P_dev(1,1,:)*mu/1e3),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('time')
    ylabel('P_d11 (kPa)')
    title('Divatoric nominal stress')

    figure
    plot(FxX,squeeze(P_dev(1,1,:)*mu/1e3),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('FxX')
    ylabel('P_d11 (kPa)')
    title('Divatoric nominal stress')
end

if pure_shear
    figure
    plot(time_series*tau1,squeeze(s_dev(1,2,:)*mu/1e3),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('time')
    ylabel('s_d12 (kPa)')
    title('Divatoric Cauchy stress')

    figure
    plot(FxY,squeeze(s_dev(1,2,:)*mu/1e3),'-o','MarkerIndices',1:markerinterval:length(time_series));
    xlabel('FxY')
    ylabel('s_d12 (kPa)')
    title('Divatoric Cauchy stress')
end