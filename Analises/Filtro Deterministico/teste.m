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


% Retirando o pedestal do conjunto de sinal de treino:

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

%% Separando o conjunto de treino (80%) e teste (20%)


ruidoDes = ruidoDes0 + ruidoDes1 + ruidoDes2 + ruidoDes3;
ruidoTes = ruidoTes0 + ruidoTes1 + ruidoTes2 + ruidoTes3;
sinalDes = sinalDes0 + sinalDes1 + sinalDes2 + sinalDes3;
sinalTes = sinalTes0 + sinalTes1 + sinalTes2 + sinalTes3;



%% scriptDayane - branqueamento -- descorrelacionando o ruido ------------------------------------
c = cov(ruidoDes);
[V,D] = eig(c); 
W = D^(-.5)*V';

%% Cálculo do pulso médio normalizado

% Filtrando os dados pra um valor de corte em MeV:

FL0C0M0 = [];

for i=1:size(sinalDes0)
    if Ma1(i,1)>1000
        FL0C0M0 = [FL0C0M0; sinalDes0(i,:)];
    end
end


FL0C1M0 = [];

for i=1:size(sinalDes1)
    if Ma1(i,2)>1000
        FL0C1M0 = [FL0C1M0; sinalDes1(i,:)];
    end
end


FL0C2M0 = [];

for i=1:size(sinalDes2)
    if Ma1(i,3)>1000
        FL0C2M0 = [FL0C2M0; sinalDes2(i,:)];
    end
end


FL0C3M0 = [];

for i=1:size(sinalDes3)
    if Ma1(i,4)>1000
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

%% Normalizando os sinais medios

sm = sm0 + sm1 + sm2 + sm3;

div = max(sm(:));
for j=1:7
    sm(j)=sm(j)/div;
end

%% PCA

NFL0M0 = [NFL0C0M0 ; NFL0C1M0 ; NFL0C2M0 ; NFL0C3M0];

[COEFF, SCORE, LATENT] = pca(sinalDes*W');

N=7;

mFC = sm*W';

mEstimacao = mFC*COEFF(:,1:N); 

%% Parametros

rRuido = (ruidoTes*W')*COEFF(:,1:N);
rSinal = (sinalTes*W')*COEFF(:,1:N);

variancia = var(ruidoDes(:,4)); 

No = variancia*2;
lambda = LATENT;

h1 = zeros(7,7);
h2 = zeros(7,7);

for i=1:N % de 1 ate o numero de pca
        h1 = h1 + ((lambda(i))./((lambda(i))+variancia))*(COEFF(:,i)*COEFF(:,i)');
        h2 = h2 + ((1./((lambda(i))+variancia)))*(COEFF(:,i)*COEFF(:,i)');
end

%% Acha a parte estocastica do sinal e do ruido

    IrRuido = zeros(size(ruidoTes,1),1); 
    IrSinal = zeros(size(sinalTes,1),1);
    for ev=1:size(ruidoTes,1) % testa qual o valor do evento usando a formula inteira do artigo
        IrRuido(ev) = (1/No)*((rRuido(ev,:)*COEFF(:,1:N)')*h1*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IrSinal(ev) = (1/No)*((rSinal(ev,:)*COEFF(:,1:N)')*h1*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

%% Acha a parte deterministica do sinal e do ruido

    IdSinal = zeros(size(sinalTes,1),1);
    IdRuido = zeros(size(ruidoTes,1),1);
    for ev=1:size(ruidoTes,1)
        IdRuido(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IdSinal(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

%% Saída do filtro

FCestRuido = IdRuido + IrRuido; % saida do filtro pra deteccao para o ruido
FCestSinal = IdSinal + IrSinal; % saida do filtro pra deteccao para o sinal

%% ROC1

pmin = min(FCestSinal);
pmax = max(FCestSinal);
pontos = 4000;

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

pmin = min(IdSinal);
pmax = max(IdSinal);
pontos = 4000;

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

%% ROC3

pmin = min(IrSinal);
pmax = max(IrSinal);
pontos = 4000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD = zeros(pontos,1);
FA = zeros(pontos,1);
pd = 0;
fa = 0;

for i=1:pontos  % patamar variando em 2000 pontos
    for j=1:size(IrRuido,1)
        if IrSinal(j) > patamar
            pd = pd + 1;
        end 
        if IrRuido(j) > patamar
            fa = fa + 1;
        end
    end
    pd = pd*100/size(IrSinal,1);
    fa = fa*100/size(IrRuido,1);
    PD(i) = pd; % preenchendo o vetor
    FA(i) = fa;
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

hold on
plot(FA, PD, '-mx');


grid
title('Análise de Eficiência - pca(sinalDes) MODULO')
legend('Filtro Estocástico (completo)', 'Componente Determinística', 'Filtro Casado', 'Componente Estocástica');
xlabel('% FA')
ylabel('% PD')



