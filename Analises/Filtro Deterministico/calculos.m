clear all; close all; clc;

load('dados\Lado-Canal.mat');
load('workspace\energia_MeV.mat');
load('noise.txt');
load('dados Lado-Canal-Modulo\L0C0M0.mat');

ruidoDes = noise(1:40308,:);
ruidoTes = noise(40309:end,:);

sinalDes = L0C0M0(1:40308,:);
sinalTes = L0C0M0(40309:end,:);

%% pedestal
pedL0C0M0 =  0;
for i=1:50385
    pedL0C0M0 = pedL0C0M0 + noise(i,1);
end

pedL0C0M0 = pedL0C0M0/50385

%% retirando o pedestal
for i=1:size(sinalDes)
    for j=1:7
        sinalDes(i,j) = sinalDes(i,j) - pedL0C0M0;
    end
end

%% Filtrando os dados pra um valor de corte em MeV

FL0C0M0 = [];

for i=1:size(sinalDes)
    if Ma1(i,1)>500
        FL0C0M0 = [FL0C0M0; sinalDes(i,:)];
    end
end

%% Normalizando os dados

NFL0C0M0 = FL0C0M0(:,:);
for i=1:size(FL0C0M0,1)
    div = max(FL0C0M0(i,:));
for j=1:7
    NFL0C0M0(i,j)=FL0C0M0(i,j)/div;
end
end

%% Pulso medio normalizado
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

medio = [a1, a2, a3, a4, a5, a6, a7];


%% Plot
figure
plot(1:7,NFL0C0M0(:,:))
title('Amostras normalizadas Lado 0 Canal 0 Modulo 0')
%axis([1 7 -5 5])
grid on

figure
plot(medio)
title('Pulso medio')
%axis([1 7 -5 5])
grid on

%% PCA

[COEFF0, SCORE0, LATENT0] = pca(NFL0C0M0);

mEstimacao = medio*COEFF0;

figure
hist(COEFF0)
title('COEFF Lado 0 Canal 0 Modulo 0')
grid

figure
plot(LATENT0,'-x')
title('LATENT Lado 0 Canal 0 Modulo 0')
grid

figure
plot(mEstimacao)
title('mEstimacao = pulsoMedio*COEFF')
grid

%% Parametros

rRuido = ruidoTes*COEFF0;
rSinal = sinalTes*COEFF0; %retirei o pedestal e filtrei pra >500MeV, nao normalizei
%rSinalNorm = NFL0C0M0 *COEFF0;

variancia = var(ruidoDes(:,4));

No = variancia*2;

lambda = LATENT0;

h1 = zeros(7,7); % vai ser a parte constante na formula de IR
h2 = zeros(7,7); % vai ser a parte constante na formula de ID

    for i=1:7 % de 1 ate o numero de pca
        h1 = h1 + ((lambda(i))./((lambda(i))+variancia))*(COEFF0(:,i)*COEFF0(:,i)');
        h2 = h2 + ((1./((lambda(i))+variancia)))*(COEFF0(:,i)*COEFF0(:,i)');
    end

figure
plot(rRuido')
title('rRuido')
grid

figure
plot(rSinal')
title('rSinal')
grid

figure
plot(h1)
title('h1')
grid

figure
plot(h2)
title('h2')
grid

%% acha a parte deterministica -----------------------
N=7;
    IdSinal = zeros(size(sinalTes,1),1);
    IdRuido = zeros(size(ruidoTes,1),1);
    for ev=1:size(ruidoTes,1)
        IdRuido(ev) = ((mEstimacao*COEFF0(:,1:N)')*h2*(rRuido(ev,:)*COEFF0(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IdSinal(ev) = ((mEstimacao*COEFF0(:,1:N)')*h2*(rSinal(ev,:)*COEFF0(:,1:N)')');
    end

    %%
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
for i=1:2000
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

stop = 1
%% Criando matrizes separando por lado, canal e modulo
L0C0M13 = [];
L0C1M13 = [];
L0C2M13 = [];
L0C3M13 = [];

m=13;

for i=m:64:3224640
   L0C0M13 = [L0C0M13; L0C0(i,4:10)];
end

for i=m:64:3224640
   L0C1M13 = [L0C1M13; L0C1(i,4:10)];
end

for i=m:64:3224640
   L0C2M13 = [L0C2M13; L0C2(i,4:10)];
end

for i=m:64:3224640
   L0C3M13 = [L0C3M13; L0C3(i,4:10)];
end

%% Filtrando os dados pra um valor de corte em MeV

FL0C0M13 = [];
FL0C1M13 = [];
FL0C2M13 = [];
FL0C3M13 = [];

for i=1:50385
    if Ma12(i,1)>150
        FL0C0M13 = [FL0C0M13; L0C0M13(i,:)];
    end
end

for i=1:50385
    if Ma12(i,2)>150
        FL0C1M13 = [FL0C1M13; L0C1M13(i,:)];
    end
end

for i=1:50385
    if Ma12(i,3)>150
        FL0C2M13 = [FL0C2M13; L0C2M13(i,:)];
    end
end

for i=1:50385
    if Ma12(i,4)>150
        FL0C3M13 = [FL0C3M13; L0C3M13(i,:)];
    end
end

%% Normalizando os dados

NFL0C0M13 = FL0C0M13(:,:);
for i=1:11901
    div = norm(NFL0C0M13(i,:));
for j=1:7
    NFL0C0M13(i,j)=NFL0C0M13(i,j)/div;
end
end

NFL0C1M13 = FL0C1M13(:,:);
for i=1:12203
    div = norm(NFL0C1M13(i,:));
for j=1:7
    NFL0C1M13(i,j)=NFL0C1M13(i,j)/div;
end
end

NFL0C2M13 = FL0C2M13(:,:);
for i=1:49606
    div = norm(NFL0C2M13(i,:));
for j=1:7
    NFL0C2M13(i,j)=NFL0C2M13(i,j)/div;
end
end

NFL0C3M13 = FL0C3M13(:,:);
for i=1:49777
    div = norm(NFL0C3M13(i,:));
for j=1:7
    NFL0C3M13(i,j)=NFL0C3M13(i,j)/div;
end
end


