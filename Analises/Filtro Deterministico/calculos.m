load('dados\Lado-Canal.mat');
load('workspace\energia_MeV.mat');
load('noise.txt');
load('dados Lado-Canal-Modulo\L0C0M0.mat');

%% pedestal
pedL0C0M0 =  0;
for i=1:50385
    pedL0C0M0 = pedL0C0M0 + noise(i,1);
end

pedL0C0M0 = pedL0C0M0/50385

%% retirando o pedestal
for i=1:50385
    for j=1:7
        L0C0M0(i,j) = L0C0M0(i,j) - pedL0C0M0;
    end
end

%% Filtrando os dados pra um valor de corte em MeV

FL0C0M0 = [];

for i=1:50385
    if Ma1(i,1)>500
        FL0C0M0 = [FL0C0M0; L0C0M0(i,:)];
    end
end

%% Normalizando os dados

NFL0C0M0 = FL0C0M0(:,:);
for i=1:8806%23344
    div = max(FL0C0M0(i,:));
for j=1:7
    NFL0C0M0(i,j)=FL0C0M0(i,j)/div;
end
end

%% Pulso medio
a1 = 0;
a2 = 0;
a3 = 0;
a4 = 0;
a5 = 0;
a6 = 0;
a7 = 0;

for i=1:8806%23344 %muda o denominador tbm
    a1 = a1 + NFL0C0M0(i,1);
    a2 = a2 + NFL0C0M0(i,2);
    a3 = a3 + NFL0C0M0(i,3);
    a4 = a4 + NFL0C0M0(i,4);
    a5 = a5 + NFL0C0M0(i,5);
    a6 = a6 + NFL0C0M0(i,6);
    a7 = a7 + NFL0C0M0(i,7);
end

a1 = a1/8806;
a2 = a2/8806;
a3 = a3/8806;
a4 = a4/8806;
a5 = a5/8806;
a6 = a6/8806;
a7 = a7/8806;

medio = [a1, a2, a3, a4, a5, a6, a7];

%% Plot
figure
plot(1:7,NFL0C0M0(:,:))
title('Lado 0 Canal 0 Modulo 0')
%axis([1 7 -5 5])
grid on

figure
plot(medio)
title('Pulso medio')
%axis([1 7 -5 5])
grid on

%% PCA

[COEFF0, SCORE0, LATENT0] = pca(FL0C0M0);

figure
hist(COEFF0)
title('COEFF Lado 0 Canal 0 Modulo 0')
grid

figure
plot(LATENT0,'-x')
title('LATENT Lado 0 Canal 0 Modulo 0')
grid

figure
plot(medio*COEFF0)
title('mEstimacao = pulsoMedio*COEFF')
grid

stop=1;

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


