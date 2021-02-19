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

%% Calculando o pedestal
ped0 =  0;
for i=1:size(noiseL0C0M0,1)
    ped0 = ped0 + noiseL0C0M0(i,1);
end
ped0 = ped0/size(noiseL0C0M0,1);


ped1 =  0;
for i=1:size(noiseL0C1M0,1)
    ped1 = ped1 + noiseL0C1M0(i,1);
end
ped1 = ped1/size(noiseL0C1M0,1);


ped2 =  0;
for i=1:size(noiseL0C2M0,1)
    ped2 = ped2 + noiseL0C2M0(i,1);
end
ped2 = ped2/size(noiseL0C2M0,1);


ped3 =  0;
for i=1:size(noiseL0C3M0,1)
    ped3 = ped3 + noiseL0C3M0(i,1);
end
ped3 = ped3/size(noiseL0C3M0,1);

%% Separando o conjunto de treino (80%) e teste (20%)

ruido0 = noiseL0C0M0;  % no scriptDayane é um pedestal mais uma distribuição gaussiana -------------
ruido1 = noiseL0C1M0;
ruido2 = noiseL0C2M0;
ruido3 = noiseL0C3M0;


ruidoDes0 = ruido0(1:40308,:);
ruidoTes0 = ruido0(40309:end,:);
sinalDes0 = L0C0M0(1:40308,:);
sinalTes0 = L0C0M0(40309:end,:);  % no scriptDayane é %ruido + distribuiçao uniforme + jitter ------


ruidoDes1 = ruido1(1:40308,:);
ruidoTes1 = ruido1(40309:end,:);
sinalDes1 = L0C1M0(1:40308,:);
sinalTes1 = L0C1M0(40309:end,:);


ruidoDes2 = ruido2(1:40308,:);
ruidoTes2 = ruido2(40309:end,:);
sinalDes2 = L0C2M0(1:40308,:);
sinalTes2 = L0C2M0(40309:end,:);


ruidoDes3 = ruido3(1:40308,:);
ruidoTes3 = ruido3(40309:end,:);
sinalDes3 = L0C3M0(1:40308,:);
sinalTes3 = L0C3M0(40309:end,:);

% %% scriptDayane - branqueamento -- descorrelacionando o ruido ------------------------------------
c = cov(ruidoDes0);
[V,D] = eig(c); 
W0 = D^(-.5)*V';

c = cov(ruidoDes0);
[V,D] = eig(c); 
W1 = D^(-.5)*V';

c = cov(ruidoDes0);
[V,D] = eig(c); 
W2 = D^(-.5)*V';

c = cov(ruidoDes0);
[V,D] = eig(c); 
W3 = D^(-.5)*V';


%% Cálculo do pulso médio normalizado

% Retirando o pedestal:

for i=1:size(sinalDes0)
    for j=1:7
        sinalDes0(i,j) = sinalDes0(i,j) - ped0;
    end
end


for i=1:size(sinalDes1)
    for j=1:7
        sinalDes1(i,j) = sinalDes1(i,j) - ped1;
    end
end


for i=1:size(sinalDes2)
    for j=1:7
        sinalDes2(i,j) = sinalDes2(i,j) - ped2;
    end
end


for i=1:size(sinalDes3)
    for j=1:7
        sinalDes3(i,j) = sinalDes3(i,j) - ped3;
    end
end


for i=1:size(ruidoDes0)
    for j=1:7
        ruidoDes0(i,j) = ruidoDes0(i,j) - ped0;
    end
end


for i=1:size(ruidoDes1)
    for j=1:7
        ruidoDes1(i,j) = ruidoDes1(i,j) - ped1;
    end
end


for i=1:size(ruidoDes2)
    for j=1:7
        ruidoDes2(i,j) = ruidoDes2(i,j) - ped2;
    end
end


for i=1:size(ruidoDes3)
    for j=1:7
        ruidoDes3(i,j) = ruidoDes3(i,j) - ped3;
    end
end

% % Retirando o pedestal do conjunto de sinal de teste:

for i=1:size(sinalTes0)
    for j=1:7
        sinalTes0(i,j) = sinalTes0(i,j) - ped0;
    end
end


for i=1:size(sinalTes1)
    for j=1:7
        sinalTes1(i,j) = sinalTes1(i,j) - ped1;
    end
end


for i=1:size(sinalTes2)
    for j=1:7
        sinalTes2(i,j) = sinalTes2(i,j) - ped2;
    end
end


for i=1:size(sinalTes3)
    for j=1:7
        sinalTes3(i,j) = sinalTes3(i,j) - ped3;
    end
end

for i=1:size(ruidoTes0)
    for j=1:7
        ruidoTes0(i,j) = ruidoTes0(i,j) - ped0;
    end
end


for i=1:size(ruidoTes1)
    for j=1:7
        ruidoTes1(i,j) = ruidoTes1(i,j) - ped1;
    end
end


for i=1:size(ruidoTes2)
    for j=1:7
        ruidoTes2(i,j) = ruidoTes2(i,j) - ped2;
    end
end


for i=1:size(ruidoTes3)
    for j=1:7
        ruidoTes3(i,j) = ruidoTes3(i,j) - ped3;
    end
end

% Filtrando os dados pra um valor de corte em MeV:

FL0C0M0 = [];

for i=1:size(sinalDes0)
    if Ma1(i,1)>500
        FL0C0M0 = [FL0C0M0; sinalDes0(i,:)];
    end
end


FL0C1M0 = [];

for i=1:size(sinalDes1)
    if Ma1(i,2)>500
        FL0C1M0 = [FL0C1M0; sinalDes1(i,:)];
    end
end


FL0C2M0 = [];

for i=1:size(sinalDes2)
    if Ma1(i,3)>500
        FL0C2M0 = [FL0C2M0; sinalDes2(i,:)];
    end
end


FL0C3M0 = [];

for i=1:size(sinalDes3)
    if Ma1(i,4)>500
        FL0C3M0 = [FL0C3M0; sinalDes3(i,:)];
    end
end

% Normalizando os dados:

NFL0C0M0 = FL0C0M0(:,:);
for i=1:size(FL0C0M0,1)
    div = max(FL0C0M0(i,:));
for j=1:7
    NFL0C0M0(i,j)=FL0C0M0(i,j)/div;
end
end


NFL0C1M0 = FL0C1M0(:,:);
for i=1:size(FL0C1M0,1)
    div = max(FL0C1M0(i,:));
for j=1:7
    NFL0C1M0(i,j)=FL0C1M0(i,j)/div;
end
end


NFL0C2M0 = FL0C2M0(:,:);
for i=1:size(FL0C2M0,1)
    div = max(FL0C2M0(i,:));
for j=1:7
    NFL0C2M0(i,j)=FL0C2M0(i,j)/div;
end
end


NFL0C3M0 = FL0C3M0(:,:);
for i=1:size(FL0C3M0,1)
    div = max(FL0C3M0(i,:));
for j=1:7
    NFL0C3M0(i,j)=FL0C3M0(i,j)/div;
end
end

% Pulso medio normalizado

a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(FL0C0M0,1)
    a1 = a1 + NFL0C0M0(i,1);
    a2 = a2 + NFL0C0M0(i,2);
    a3 = a3 + NFL0C0M0(i,3);
    a4 = a4 + NFL0C0M0(i,4);
    a5 = a5 + NFL0C0M0(i,5);
    a6 = a6 + NFL0C0M0(i,6);
    a7 = a7 + NFL0C0M0(i,7);
end
a1 = a1/size(FL0C0M0,1);
a2 = a2/size(FL0C0M0,1);
a3 = a3/size(FL0C0M0,1);
a4 = a4/size(FL0C0M0,1);
a5 = a5/size(FL0C0M0,1);
a6 = a6/size(FL0C0M0,1);
a7 = a7/size(FL0C0M0,1);
sm0 = [a1, a2, a3, a4, a5, a6, a7];


a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(FL0C1M0,1)
    a1 = a1 + NFL0C1M0(i,1);
    a2 = a2 + NFL0C1M0(i,2);
    a3 = a3 + NFL0C1M0(i,3);
    a4 = a4 + NFL0C1M0(i,4);
    a5 = a5 + NFL0C1M0(i,5);
    a6 = a6 + NFL0C1M0(i,6);
    a7 = a7 + NFL0C1M0(i,7);
end
a1 = a1/size(FL0C1M0,1);
a2 = a2/size(FL0C1M0,1);
a3 = a3/size(FL0C1M0,1);
a4 = a4/size(FL0C1M0,1);
a5 = a5/size(FL0C1M0,1);
a6 = a6/size(FL0C1M0,1);
a7 = a7/size(FL0C1M0,1);
sm1 = [a1, a2, a3, a4, a5, a6, a7];


a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(FL0C2M0,1)
    a1 = a1 + NFL0C2M0(i,1);
    a2 = a2 + NFL0C2M0(i,2);
    a3 = a3 + NFL0C2M0(i,3);
    a4 = a4 + NFL0C2M0(i,4);
    a5 = a5 + NFL0C2M0(i,5);
    a6 = a6 + NFL0C2M0(i,6);
    a7 = a7 + NFL0C2M0(i,7);
end
a1 = a1/size(FL0C2M0,1);
a2 = a2/size(FL0C2M0,1);
a3 = a3/size(FL0C2M0,1);
a4 = a4/size(FL0C2M0,1);
a5 = a5/size(FL0C2M0,1);
a6 = a6/size(FL0C2M0,1);
a7 = a7/size(FL0C2M0,1);
sm2 = [a1, a2, a3, a4, a5, a6, a7];


a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;
for i=1:size(FL0C3M0,1)
    a1 = a1 + NFL0C3M0(i,1);
    a2 = a2 + NFL0C3M0(i,2);
    a3 = a3 + NFL0C3M0(i,3);
    a4 = a4 + NFL0C3M0(i,4);
    a5 = a5 + NFL0C3M0(i,5);
    a6 = a6 + NFL0C3M0(i,6);
    a7 = a7 + NFL0C3M0(i,7);
end
a1 = a1/size(FL0C3M0,1);
a2 = a2/size(FL0C3M0,1);
a3 = a3/size(FL0C3M0,1);
a4 = a4/size(FL0C3M0,1);
a5 = a5/size(FL0C3M0,1);
a6 = a6/size(FL0C3M0,1);
a7 = a7/size(FL0C3M0,1);
sm3 = [a1, a2, a3, a4, a5, a6, a7];

%% Plot
% figure
% subplot(2,2,1);
% plot(1:7,NFL0C0M0(:,:))
% title('Amostras normalizadas Canal 0 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,2);
% plot(1:7,NFL0C1M0(:,:))
% title('Amostras normalizadas Canal 1 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,3);
% plot(1:7,NFL0C2M0(:,:))
% title('Amostras normalizadas Canal 2 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,4);
% plot(1:7,NFL0C3M0(:,:))
% title('Amostras normalizadas Canal 3 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% figure
% subplot(2,2,1);
% plot(sm0)
% title('Pulso medio normalizado Canal 0 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,2);
% plot(sm1)
% title('Pulso medio normalizado Canal 1 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,3);
% plot(sm2)
% title('Pulso medio normalizado Canal 2 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,4);
% plot(sm3)
% title('Pulso medio normalizado Canal 3 Modulo 0')
% %axis([1 7 -5 5])
% grid on

%% PCA

[COEFF0, SCORE0, LATENT0] = pca(NFL0C0M0*W0');  % scriptDayane usou o pegaPulseJitter()*W'--------------
[COEFF1, SCORE1, LATENT1] = pca(NFL0C1M0*W1');
[COEFF2, SCORE2, LATENT2] = pca(NFL0C2M0*W2');
[COEFF3, SCORE3, LATENT3] = pca(NFL0C3M0*W3');


N=7;

mFC0 = sm0*W0';
mFC1 = sm1*W1';
mFC2 = sm2*W2';
mFC3 = sm3*W3';

mEstimacao0 = mFC0*COEFF0(:,1:N); % scriptDayane usou sm*W' branqueado ------------------------------
mEstimacao1 = mFC1*COEFF1(:,1:N);
mEstimacao2 = mFC2*COEFF2(:,1:N);
mEstimacao3 = mFC3*COEFF3(:,1:N);

%% Plot
% figure
% subplot(2,2,1);
% plot(LATENT0,'-x')
% title('LATENT Canal 0 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,2);
% plot(LATENT1,'-x')
% title('LATENT Canal 1 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,3);
% plot(LATENT2,'-x')
% title('LATENT Canal 2 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% subplot(2,2,4);
% plot(LATENT3,'-x')
% title('LATENT Canal 3 Modulo 0')
% %axis([1 7 -5 5])
% grid on

%% Parametros

rRuido0 = (ruidoTes0*W0')*COEFF0(:,1:N);
rSinal0 = (sinalTes0*W0')*COEFF0(:,1:N); % SmpMuon puro. no scriptDayane subtrai o pedestal ------------ 

rRuido1 = (ruidoTes1*W1')*COEFF1(:,1:N);
rSinal1 = (sinalTes1*W1')*COEFF1(:,1:N);

rRuido2 = (ruidoTes2*W2')*COEFF2(:,1:N);
rSinal2 = (sinalTes2*W2')*COEFF2(:,1:N);

rRuido3 = (ruidoTes3*W3')*COEFF3(:,1:N);
rSinal3 = (sinalTes3*W3')*COEFF3(:,1:N);

variancia0 = var(ruidoDes0(:,4)); % usei o mesmo parametro do scriptDayane -------------------------

variancia1 = var(ruidoDes1(:,4));

variancia2 = var(ruidoDes2(:,4));

variancia3 = var(ruidoDes3(:,4));

No0 = variancia0*2;
lambda0 = LATENT0;

No1 = variancia1*2;
lambda1 = LATENT1;

No2 = variancia2*2;
lambda2 = LATENT2;

No3 = variancia3*2;
lambda3 = LATENT3;

h10 = zeros(7,7);
h20 = zeros(7,7);

h11 = zeros(7,7);
h21 = zeros(7,7);

h12 = zeros(7,7);
h22 = zeros(7,7);

h13 = zeros(7,7);
h23 = zeros(7,7);

    for i=1:7 % de 1 ate o numero de pca
        h10 = h10 + ((lambda0(i))./((lambda0(i))+variancia0))*(COEFF0(:,i)*COEFF0(:,i)');
        h20 = h20 + ((1./((lambda0(i))+variancia0)))*(COEFF0(:,i)*COEFF0(:,i)');
        
        h11 = h11 + ((lambda1(i))./((lambda1(i))+variancia1))*(COEFF1(:,i)*COEFF1(:,i)');
        h21 = h21 + ((1./((lambda1(i))+variancia1)))*(COEFF1(:,i)*COEFF1(:,i)');
        
        h12 = h12 + ((lambda2(i))./((lambda2(i))+variancia2))*(COEFF2(:,i)*COEFF2(:,i)');
        h22 = h22 + ((1./((lambda2(i))+variancia2)))*(COEFF2(:,i)*COEFF2(:,i)');
        
        h13 = h13 + ((lambda3(i))./((lambda3(i))+variancia3))*(COEFF3(:,i)*COEFF3(:,i)');
        h23 = h23 + ((1./((lambda3(i))+variancia3)))*(COEFF3(:,i)*COEFF3(:,i)');
    end

% %% Plot
% figure
% plot(rRuido')
% title('rRuido')
% grid
% 
% figure
% plot(rSinal')
% title('rSinal')
% grid
% 
% figure
% plot(h1)
% title('h1')
% grid
% 
% figure
% plot(h2)
% title('h2')
% grid

    
%% Acha a parte deterministica do sinal e do ruido

N=7;

IdSinal0 = zeros(size(sinalTes0,1),1);
IdRuido0 = zeros(size(ruidoTes0,1),1);
for ev=1:size(ruidoTes0,1)
    IdRuido0(ev) = ((mEstimacao0*COEFF0(:,1:N)')*h20*(rRuido0(ev,:)*COEFF0(:,1:N)')');
end

for ev=1:size(sinalTes0,1)
    IdSinal0(ev) = ((mEstimacao0*COEFF0(:,1:N)')*h20*(rSinal0(ev,:)*COEFF0(:,1:N)')');
end


IdSinal1 = zeros(size(sinalTes1,1),1);
IdRuido1 = zeros(size(ruidoTes1,1),1);
for ev=1:size(ruidoTes1,1)
    IdRuido1(ev) = ((mEstimacao1*COEFF1(:,1:N)')*h21*(rRuido1(ev,:)*COEFF1(:,1:N)')');
end

for ev=1:size(sinalTes1,1)
    IdSinal1(ev) = ((mEstimacao1*COEFF1(:,1:N)')*h21*(rSinal1(ev,:)*COEFF1(:,1:N)')');
end


IdSinal2 = zeros(size(sinalTes2,1),1);
IdRuido2 = zeros(size(ruidoTes2,1),1);
for ev=1:size(ruidoTes2,1)
    IdRuido2(ev) = ((mEstimacao2*COEFF2(:,1:N)')*h22*(rRuido2(ev,:)*COEFF2(:,1:N)')');
end

for ev=1:size(sinalTes2,1)
    IdSinal2(ev) = ((mEstimacao2*COEFF2(:,1:N)')*h22*(rSinal2(ev,:)*COEFF2(:,1:N)')');
end


IdSinal3 = zeros(size(sinalTes3,1),1);
IdRuido3 = zeros(size(ruidoTes3,1),1);
for ev=1:size(ruidoTes3,1)
    IdRuido3(ev) = ((mEstimacao3*COEFF3(:,1:N)')*h23*(rRuido3(ev,:)*COEFF3(:,1:N)')');
end

for ev=1:size(sinalTes3,1)
    IdSinal3(ev) = ((mEstimacao3*COEFF3(:,1:N)')*h23*(rSinal3(ev,:)*COEFF3(:,1:N)')');
end

%% Compondo os sinais

IdRuido = IdRuido0 + IdRuido1 + IdRuido2 + IdRuido3;
IdSinal = IdSinal0 + IdSinal1 + IdSinal2 + IdSinal3;

%% Acha a parte estocastica do sinal e do ruido

    IrRuido0 = zeros(size(ruidoTes0,1),1); 
    IrSinal0 = zeros(size(sinalTes0,1),1);
    
    IrRuido1 = zeros(size(ruidoTes1,1),1); 
    IrSinal1 = zeros(size(sinalTes1,1),1);
    
    IrRuido2 = zeros(size(ruidoTes2,1),1); 
    IrSinal2 = zeros(size(sinalTes2,1),1);
    
    IrRuido3 = zeros(size(ruidoTes3,1),1); 
    IrSinal3 = zeros(size(sinalTes3,1),1);
    
    
    for ev=1:size(ruidoTes0,1)
        IrRuido0(ev) = (1/No0)*((rRuido0(ev,:)*COEFF0(:,1:N)')*h10*(rRuido0(ev,:)*COEFF0(:,1:N)')');
    end
    for ev=1:size(sinalTes0,1)
        IrSinal0(ev) = (1/No0)*((rSinal0(ev,:)*COEFF0(:,1:N)')*h10*(rSinal0(ev,:)*COEFF0(:,1:N)')');
    end
    
    
    for ev=1:size(ruidoTes1,1)
        IrRuido1(ev) = (1/No1)*((rRuido1(ev,:)*COEFF1(:,1:N)')*h11*(rRuido1(ev,:)*COEFF1(:,1:N)')');
    end
    for ev=1:size(sinalTes1,1)
        IrSinal1(ev) = (1/No1)*((rSinal1(ev,:)*COEFF1(:,1:N)')*h11*(rSinal1(ev,:)*COEFF1(:,1:N)')');
    end
    
    
    for ev=1:size(ruidoTes2,1)
        IrRuido2(ev) = (1/No2)*((rRuido2(ev,:)*COEFF2(:,1:N)')*h12*(rRuido2(ev,:)*COEFF2(:,1:N)')');
    end
    for ev=1:size(sinalTes2,1)
        IrSinal2(ev) = (1/No2)*((rSinal2(ev,:)*COEFF2(:,1:N)')*h12*(rSinal2(ev,:)*COEFF2(:,1:N)')');
    end
    
    
    for ev=1:size(ruidoTes3,1)
        IrRuido3(ev) = (1/No3)*((rRuido3(ev,:)*COEFF3(:,1:N)')*h13*(rRuido3(ev,:)*COEFF3(:,1:N)')');
    end
    for ev=1:size(sinalTes3,1)
        IrSinal3(ev) = (1/No3)*((rSinal3(ev,:)*COEFF3(:,1:N)')*h13*(rSinal3(ev,:)*COEFF3(:,1:N)')');
    end

%% Compondo os sinais

IrRuido = IrRuido0 + IrRuido1 + IrRuido2 + IrRuido3;
IrSinal = IrSinal0 + IrSinal1 + IrSinal2 + IrSinal3;

%% Saída do filtro

FCestRuido = IdRuido + IrRuido; % saida do filtro pra deteccao para o ruido
FCestSinal = IdSinal + IrSinal; % saida do filtro pra deteccao para o sinal

%% Plot
% figure
% plot(FCestRuido)
% title('FCestRuido')
% grid
% 
% figure
% plot(FCestSinal)
% title('FCestSinal')
% grid

%% ROC1

pmin = -2;
pmax = 2;
pontos = 1000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD1 = zeros(pontos,1);
FA1 = zeros(pontos,1);
pd = 0;
fa = 0;

for i=1:pontos  % patamar variando em x pontos
    for j=1:size(FCestRuido,1)
        if FCestSinal(j) > patamar
            pd = pd + 1;
        end 
        if FCestRuido(j) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(FCestSinal,1);
    fa = fa*100/size(FCestRuido,1);
    PD1(i) = pd; % preenchendo o vetor
    FA1(i) = fa;
    pd = 0;
    fa = 0;
    patamar = patamar + psoma;
end

%% ROC2

pmin = -2;
pmax = 20;
pontos = 2000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD2 = zeros(pontos,1);
FA2 = zeros(pontos,1);
pd = 0;
fa = 0;

for i=1:pontos  % patamar variando em 2000 pontos
    for j=1:size(IdRuido,1)
        if IdSinal(j) > patamar
            pd = pd + 1;
        end 
        if IdRuido(j) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(IdSinal,1);
    fa = fa*100/size(IdRuido,1);
    PD2(i) = pd; % preenchendo o vetor
    FA2(i) = fa;
    pd = 0;
    fa = 0;
    patamar = patamar + psoma;
end

%% ROC Filtro Casado


IdSinal = [];
IdRuido = [];

IdSinal = muonMa1(:,1) + muonMa1(:,2) + muonMa1(:,3) + muonMa1(:,4);
IdRuido = noiseMa1(:,1) + noiseMa1(:,2) + noiseMa1(:,3) + noiseMa1(:,4);

pmin = -1020;
pmax = 3020;
pontos = 1000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD3 = zeros(pontos,1);
FA3 = zeros(pontos,1);
pd = 0;
fa = 0;

for i=1:pontos  % patamar variando em 2000 pontos
    for j=1:size(IdRuido,1)
        if IdSinal(j) > patamar
            pd = pd + 1;
        end 
        if IdRuido(j) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(IdSinal,1);
    fa = fa*100/size(IdRuido,1);
    PD3(i) = pd; % preenchendo o vetor
    FA3(i) = fa;
    pd = 0;
    fa = 0;
    patamar = patamar + psoma;
end


figure

plot(FA1, PD1, '-b*');

hold on
plot(FA2, PD2, '-rs');

hold on
plot(FA3, PD3, '-.g', 'LineWidth', 2);


grid
title('Análise de Eficiência - ROC')
legend('Filtro Estocástico (completo)', 'Filtro Estocástico (simplificado)', 'Filtro Casado');
xlabel('% FA')
ylabel('% PD')



