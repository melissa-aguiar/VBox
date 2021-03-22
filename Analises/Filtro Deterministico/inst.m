clear all;
close all;
clc;

%% Carregando arquivos

load('workspace\energia_MeV.mat');
load('noiseL0C0M0.txt');
load('noiseL0C1M0.txt');
load('noiseL0C2M0.txt');
load('noiseL0C3M0.txt');
load('dados Lado-Canal-Modulo\L0C0M0.mat');
load('dados Lado-Canal-Modulo\L0C1M0.mat');
load('dados Lado-Canal-Modulo\L0C2M0.mat');
load('dados Lado-Canal-Modulo\L0C3M0.mat');

load('muonMa1.mat');
load('noiseMa1.mat');

%% Plot do ruido
figure
subplot(2,2,1);
plot(noiseL0C0M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Channel 0')

subplot(2,2,2);
plot(noiseL0C1M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Channel 1')

subplot(2,2,3);
plot(noiseL0C2M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Channel 2')

subplot(2,2,4);
plot(noiseL0C3M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Channel 3')

%% Plot do sinal
figure
subplot(2,2,1);
plot(L0C0M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 0')

subplot(2,2,2);
plot(L0C1M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 1')

subplot(2,2,3);
plot(L0C2M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 2')

subplot(2,2,4);
plot(L0C3M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 3')

%% Pulso medio ruido

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(noiseL0C0M0,1)
    a1 = a1 + noiseL0C0M0(i,1);
    a2 = a2 + noiseL0C0M0(i,2);
    a3 = a3 + noiseL0C0M0(i,3);
    a4 = a4 + noiseL0C0M0(i,4);
    a5 = a5 + noiseL0C0M0(i,5);
    a6 = a6 + noiseL0C0M0(i,6);
    a7 = a7 + noiseL0C0M0(i,7);
end
a1 = a1/size(noiseL0C0M0,1);
a2 = a2/size(noiseL0C0M0,1);
a3 = a3/size(noiseL0C0M0,1);
a4 = a4/size(noiseL0C0M0,1);
a5 = a5/size(noiseL0C0M0,1);
a6 = a6/size(noiseL0C0M0,1);
a7 = a7/size(noiseL0C0M0,1);
sm0 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(noiseL0C1M0,1)
    a1 = a1 + noiseL0C1M0(i,1);
    a2 = a2 + noiseL0C1M0(i,2);
    a3 = a3 + noiseL0C1M0(i,3);
    a4 = a4 + noiseL0C1M0(i,4);
    a5 = a5 + noiseL0C1M0(i,5);
    a6 = a6 + noiseL0C1M0(i,6);
    a7 = a7 + noiseL0C1M0(i,7);
end
a1 = a1/size(noiseL0C1M0,1);
a2 = a2/size(noiseL0C1M0,1);
a3 = a3/size(noiseL0C1M0,1);
a4 = a4/size(noiseL0C1M0,1);
a5 = a5/size(noiseL0C1M0,1);
a6 = a6/size(noiseL0C1M0,1);
a7 = a7/size(noiseL0C1M0,1);
sm1 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(noiseL0C2M0,1)
    a1 = a1 + noiseL0C2M0(i,1);
    a2 = a2 + noiseL0C2M0(i,2);
    a3 = a3 + noiseL0C2M0(i,3);
    a4 = a4 + noiseL0C2M0(i,4);
    a5 = a5 + noiseL0C2M0(i,5);
    a6 = a6 + noiseL0C2M0(i,6);
    a7 = a7 + noiseL0C2M0(i,7);
end
a1 = a1/size(noiseL0C2M0,1);
a2 = a2/size(noiseL0C2M0,1);
a3 = a3/size(noiseL0C2M0,1);
a4 = a4/size(noiseL0C2M0,1);
a5 = a5/size(noiseL0C2M0,1);
a6 = a6/size(noiseL0C2M0,1);
a7 = a7/size(noiseL0C2M0,1);
sm2 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(noiseL0C3M0,1)
    a1 = a1 + noiseL0C3M0(i,1);
    a2 = a2 + noiseL0C3M0(i,2);
    a3 = a3 + noiseL0C3M0(i,3);
    a4 = a4 + noiseL0C3M0(i,4);
    a5 = a5 + noiseL0C3M0(i,5);
    a6 = a6 + noiseL0C3M0(i,6);
    a7 = a7 + noiseL0C3M0(i,7);
end
a1 = a1/size(noiseL0C3M0,1);
a2 = a2/size(noiseL0C3M0,1);
a3 = a3/size(noiseL0C3M0,1);
a4 = a4/size(noiseL0C3M0,1);
a5 = a5/size(noiseL0C3M0,1);
a6 = a6/size(noiseL0C3M0,1);
a7 = a7/size(noiseL0C3M0,1);
sm3 = [a1, a2, a3, a4, a5, a6, a7];

figure
subplot(2,2,1);
plot(sm0,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Mean Values Channel 0')

subplot(2,2,2);
plot(sm1,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Mean Values Channel 1')

subplot(2,2,3);
plot(sm2,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Mean Values Channel 2')

subplot(2,2,4);
plot(sm3,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Noise Mean Values Channel 3')

%% Pulso medio sinal

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C0M0,1)
    a1 = a1 + L0C0M0(i,1);
    a2 = a2 + L0C0M0(i,2);
    a3 = a3 + L0C0M0(i,3);
    a4 = a4 + L0C0M0(i,4);
    a5 = a5 + L0C0M0(i,5);
    a6 = a6 + L0C0M0(i,6);
    a7 = a7 + L0C0M0(i,7);
end
a1 = a1/size(L0C0M0,1);
a2 = a2/size(L0C0M0,1);
a3 = a3/size(L0C0M0,1);
a4 = a4/size(L0C0M0,1);
a5 = a5/size(L0C0M0,1);
a6 = a6/size(L0C0M0,1);
a7 = a7/size(L0C0M0,1);
sm0 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C1M0,1)
    a1 = a1 + L0C1M0(i,1);
    a2 = a2 + L0C1M0(i,2);
    a3 = a3 + L0C1M0(i,3);
    a4 = a4 + L0C1M0(i,4);
    a5 = a5 + L0C1M0(i,5);
    a6 = a6 + L0C1M0(i,6);
    a7 = a7 + L0C1M0(i,7);
end
a1 = a1/size(L0C1M0,1);
a2 = a2/size(L0C1M0,1);
a3 = a3/size(L0C1M0,1);
a4 = a4/size(L0C1M0,1);
a5 = a5/size(L0C1M0,1);
a6 = a6/size(L0C1M0,1);
a7 = a7/size(L0C1M0,1);
sm1 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C2M0,1)
    a1 = a1 + L0C2M0(i,1);
    a2 = a2 + L0C2M0(i,2);
    a3 = a3 + L0C2M0(i,3);
    a4 = a4 + L0C2M0(i,4);
    a5 = a5 + L0C2M0(i,5);
    a6 = a6 + L0C2M0(i,6);
    a7 = a7 + L0C2M0(i,7);
end
a1 = a1/size(L0C2M0,1);
a2 = a2/size(L0C2M0,1);
a3 = a3/size(L0C2M0,1);
a4 = a4/size(L0C2M0,1);
a5 = a5/size(L0C2M0,1);
a6 = a6/size(L0C2M0,1);
a7 = a7/size(L0C2M0,1);
sm2 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C3M0,1)
    a1 = a1 + L0C3M0(i,1);
    a2 = a2 + L0C3M0(i,2);
    a3 = a3 + L0C3M0(i,3);
    a4 = a4 + L0C3M0(i,4);
    a5 = a5 + L0C3M0(i,5);
    a6 = a6 + L0C3M0(i,6);
    a7 = a7 + L0C3M0(i,7);
end
a1 = a1/size(L0C3M0,1);
a2 = a2/size(L0C3M0,1);
a3 = a3/size(L0C3M0,1);
a4 = a4/size(L0C3M0,1);
a5 = a5/size(L0C3M0,1);
a6 = a6/size(L0C3M0,1);
a7 = a7/size(L0C3M0,1);
sm3 = [a1, a2, a3, a4, a5, a6, a7];

figure
subplot(2,2,1);
plot(sm0,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 0')

subplot(2,2,2);
plot(sm1,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 1')

subplot(2,2,3);
plot(sm2,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 2')

subplot(2,2,4);
plot(sm3,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 3')

%% Plot do sinal sem pedestal
figure
subplot(2,2,1);
plot(L0C0M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 0')

subplot(2,2,2);
plot(L0C1M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 1')

subplot(2,2,3);
plot(L0C2M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 2')

subplot(2,2,4);
plot(L0C3M0(1:50385,:)')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Channel 3')

%% Pulso médio do sinal sem pedestal

% Calculando o pedestal
ped0 =  0;
for i=1:size(noiseL0C0M0,1)
    ped0 = ped0 + noiseL0C0M0(i,4);
end
ped0 = ped0/size(noiseL0C0M0,1);


ped1 =  0;
for i=1:size(noiseL0C1M0,1)
    ped1 = ped1 + noiseL0C1M0(i,4);
end
ped1 = ped1/size(noiseL0C1M0,1);


ped2 =  0;
for i=1:size(noiseL0C2M0,1)
    ped2 = ped2 + noiseL0C2M0(i,4);
end
ped2 = ped2/size(noiseL0C2M0,1);


ped3 =  0;
for i=1:size(noiseL0C3M0,1)
    ped3 = ped3 + noiseL0C3M0(i,4);
end
ped3 = ped3/size(noiseL0C3M0,1);

% Retirando o pedestal do conjunto de sinal:

for i=1:size(L0C0M0)
    for j=1:7
        L0C0M0(i,j) = L0C0M0(i,j) - ped0;
    end
end

for i=1:size(L0C1M0)
    for j=1:7
        L0C1M0(i,j) = L0C1M0(i,j) - ped0;
    end
end

for i=1:size(L0C2M0)
    for j=1:7
        L0C2M0(i,j) = L0C2M0(i,j) - ped0;
    end
end

for i=1:size(L0C3M0)
    for j=1:7
        L0C3M0(i,j) = L0C3M0(i,j) - ped0;
    end
end

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C0M0,1)
    a1 = a1 + L0C0M0(i,1);
    a2 = a2 + L0C0M0(i,2);
    a3 = a3 + L0C0M0(i,3);
    a4 = a4 + L0C0M0(i,4);
    a5 = a5 + L0C0M0(i,5);
    a6 = a6 + L0C0M0(i,6);
    a7 = a7 + L0C0M0(i,7);
end
a1 = a1/size(L0C0M0,1);
a2 = a2/size(L0C0M0,1);
a3 = a3/size(L0C0M0,1);
a4 = a4/size(L0C0M0,1);
a5 = a5/size(L0C0M0,1);
a6 = a6/size(L0C0M0,1);
a7 = a7/size(L0C0M0,1);
sm0 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C1M0,1)
    a1 = a1 + L0C1M0(i,1);
    a2 = a2 + L0C1M0(i,2);
    a3 = a3 + L0C1M0(i,3);
    a4 = a4 + L0C1M0(i,4);
    a5 = a5 + L0C1M0(i,5);
    a6 = a6 + L0C1M0(i,6);
    a7 = a7 + L0C1M0(i,7);
end
a1 = a1/size(L0C1M0,1);
a2 = a2/size(L0C1M0,1);
a3 = a3/size(L0C1M0,1);
a4 = a4/size(L0C1M0,1);
a5 = a5/size(L0C1M0,1);
a6 = a6/size(L0C1M0,1);
a7 = a7/size(L0C1M0,1);
sm1 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C2M0,1)
    a1 = a1 + L0C2M0(i,1);
    a2 = a2 + L0C2M0(i,2);
    a3 = a3 + L0C2M0(i,3);
    a4 = a4 + L0C2M0(i,4);
    a5 = a5 + L0C2M0(i,5);
    a6 = a6 + L0C2M0(i,6);
    a7 = a7 + L0C2M0(i,7);
end
a1 = a1/size(L0C2M0,1);
a2 = a2/size(L0C2M0,1);
a3 = a3/size(L0C2M0,1);
a4 = a4/size(L0C2M0,1);
a5 = a5/size(L0C2M0,1);
a6 = a6/size(L0C2M0,1);
a7 = a7/size(L0C2M0,1);
sm2 = [a1, a2, a3, a4, a5, a6, a7];

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(L0C3M0,1)
    a1 = a1 + L0C3M0(i,1);
    a2 = a2 + L0C3M0(i,2);
    a3 = a3 + L0C3M0(i,3);
    a4 = a4 + L0C3M0(i,4);
    a5 = a5 + L0C3M0(i,5);
    a6 = a6 + L0C3M0(i,6);
    a7 = a7 + L0C3M0(i,7);
end
a1 = a1/size(L0C3M0,1);
a2 = a2/size(L0C3M0,1);
a3 = a3/size(L0C3M0,1);
a4 = a4/size(L0C3M0,1);
a5 = a5/size(L0C3M0,1);
a6 = a6/size(L0C3M0,1);
a7 = a7/size(L0C3M0,1);
sm3 = [a1, a2, a3, a4, a5, a6, a7];


figure
subplot(2,2,1);
plot(sm0,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 0')

subplot(2,2,2);
plot(sm1,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 1')

subplot(2,2,3);
plot(sm2,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 2')

subplot(2,2,4);
plot(sm3,'-x')
grid
xlabel('Sample')
ylabel('Amplitude')
title('Signal Mean Values Channel 3')

%% ROC Filtro Casado


FcSinal = muonMa1(:,1) + muonMa1(:,2) + muonMa1(:,3) + muonMa1(:,4);
FcRuido = noiseMa1(:,1) + noiseMa1(:,2) + noiseMa1(:,3) + noiseMa1(:,4);

pmin = min(FcSinal);
pmax = max(FcSinal);
pontos = 1000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD3 = zeros(pontos,1);
FA3 = zeros(pontos,1);
pd = 0;
fa = 0;

for i=1:pontos  % patamar variando em 2000 pontos
    for j=1:size(FcRuido,1)
        if FcSinal(j) > patamar
            pd = pd + 1;
        end 
        if FcRuido(j) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(FcSinal,1);
    fa = fa*100/size(FcRuido,1);
    PD3(i) = pd; % preenchendo o vetor
    FA3(i) = fa;
    pd = 0;
    fa = 0;
    patamar = patamar + psoma;
end
