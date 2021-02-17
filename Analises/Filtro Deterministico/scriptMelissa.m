clear all;
close all;
clc;

%% Carregando arquivos

load('dados\Lado-Canal.mat');
load('workspace\energia_MeV.mat');
load('noise.txt');
load('dados Lado-Canal-Modulo\L0C0M0.mat');

%% Calculando o pedestal
ped =  0;
for i=1:size(noise,1)
    ped = ped + noise(i,1);
end

ped = ped/size(noise,1);

%% Separando o conjunto de treino (80%) e teste (20%)

ruido = noise;  % é o SmpNoise[0][0][0] ------------ no scriptDayane é um pedestal mais uma distribuição gaussiana -----

ruidoDes = ruido(1:40308,:);
ruidoTes = ruido(40309:end,:);

sinalDes = L0C0M0(1:40308,:);
sinalTes = L0C0M0(40309:end,:);  % no scriptDayane é %ruido + distribuiçao uniforme + jitter ---------------------------

% %% scriptDayane - branqueamento -- descorrelacionando o ruido --------------------------------------------------------
% c = cov(ruidoDes); %ruidoDes é descorrelacionado
% [V,D] = eig(c); % autovalores e autovetores
% W = D^(-.5)*V'; % branqueamento
% %ruido branqueado

%% Cálculo do pulso médio normalizado

% Retirando o pedestal:

for i=1:size(sinalDes)
    for j=1:7
        sinalDes(i,j) = sinalDes(i,j) - ped;
    end
end

% Filtrando os dados pra um valor de corte em MeV:

FL0C0M0 = [];

for i=1:size(sinalDes)
    if Ma1(i,1)>500
        FL0C0M0 = [FL0C0M0; sinalDes(i,:)];
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

sm = [a1, a2, a3, a4, a5, a6, a7];

% %% Plot
% figure
% plot(1:7,NFL0C0M0(:,:))
% title('Amostras normalizadas Lado 0 Canal 0 Modulo 0')
% %axis([1 7 -5 5])
% grid on
% 
% figure
% plot(medio)
% title('Pulso medio')
% %axis([1 7 -5 5])
% grid on

%% PCA

[COEFF0, SCORE0, LATENT0] = pca(NFL0C0M0);  % scriptDayane usou o pegaPulseJitter()*W'----------------------------------

N=7;
mEstimacao = sm*COEFF0(:,1:N); % scriptDayane usou sm*W' branqueado ----------------------------------------------------

% %% Plot
% figure
% hist(COEFF0)
% title('COEFF Lado 0 Canal 0 Modulo 0')
% grid
% 
% figure
% plot(LATENT0,'-x')
% title('LATENT Lado 0 Canal 0 Modulo 0')
% grid
% 
% figure
% plot(mEstimacao)
% title('mEstimacao = pulsoMedio*COEFF')
% grid

%% Parametros

rRuido = (ruidoTes)*COEFF0(:,1:N);
rSinal = (sinalTes)*COEFF0(:,1:N); % SmpMuon puro.. no scriptDayane o pedestal é retirado do sinal ---------------------
%rSinalNorm = NFL0C0M0 *COEFF0;

variancia = var(ruidoDes(:,4)); % usei o mesmo parametro do scriptDayane -----------------------------------------------

No = variancia*2;

lambda = LATENT0;

h1 = zeros(7,7); % vai ser a parte constante na formula de IR
h2 = zeros(7,7); % vai ser a parte constante na formula de ID

    for i=1:7 % de 1 ate o numero de pca
        h1 = h1 + ((lambda(i))./((lambda(i))+variancia))*(COEFF0(:,i)*COEFF0(:,i)');
        h2 = h2 + ((1./((lambda(i))+variancia)))*(COEFF0(:,i)*COEFF0(:,i)');
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
IdSinal = zeros(size(sinalTes,1),1);
IdRuido = zeros(size(ruidoTes,1),1);
for ev=1:size(ruidoTes,1)
    IdRuido(ev) = ((mEstimacao*COEFF0(:,1:N)')*h2*(rRuido(ev,:)*COEFF0(:,1:N)')');
end

for ev=1:size(sinalTes,1)
    IdSinal(ev) = ((mEstimacao*COEFF0(:,1:N)')*h2*(rSinal(ev,:)*COEFF0(:,1:N)')');
end

%% Plot
figure
plot(IdRuido)
title('IdRuido')
grid

figure
plot(IdSinal)
title('IdSinal')
grid

%% ROC

patamar = 0;
PD = [];
FA = [];
pd = 0;
fa = 0;

for i=1:2000  % patamar variando em 2000 pontos
    for j=1:size(IdRuido,1)
        if IdSinal(j,1) > patamar
            pd = pd + 1;
        end
        if IdRuido(j,1) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(IdRuido,1);
    fa = fa*100/size(IdRuido,1);
    PD = [PD; pd];
    FA = [FA; fa];
    pd = 0;
    fa = 0;
    patamar = patamar + 0.01;
end

figure
plot(FA, PD, '-x')
grid
title('ROC')
xlabel('% FA')
ylabel('% PD')

