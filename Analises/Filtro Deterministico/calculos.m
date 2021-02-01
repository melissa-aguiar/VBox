%% Plotando os pulsos Lado0-Modulo0-Canais0123

figure
plot(1:7,NFL0C0M0(:,:))
title('Lado 0 Canal 0 Modulo 0')
grid on

figure
plot(1:7,NFL0C1M0(:,:))
title('Lado 0 Canal 1 Modulo 0')
grid on

figure
plot(1:7,NFL0C2M0(:,:))
title('Lado 0 Canal 2 Modulo 0')
grid on

figure
plot(1:7,NFL0C3M0(:,:))
title('Lado 0 Canal 3 Modulo 0')
grid on

stop=0;


clear all
close all
clc

tic

%% Criando matrizes separando por lado, canal e modulo
load('dados\Lado-Canal.mat');
L0C0M0 = [];
L0C1M0 = [];
L0C2M0 = [];
L0C3M0 = [];

m=1;

for i=m:64:3224640
   L0C0M0 = [L0C0M0; L0C0(i,4:10)];
end

for i=m:64:3224640
   L0C1M0 = [L0C1M0; L0C1(i,4:10)];
end

for i=m:64:3224640
   L0C2M0 = [L0C2M0; L0C2(i,4:10)];
end

for i=m:64:3224640
   L0C3M0 = [L0C3M0; L0C3(i,4:10)];
end

%% Filtrando os dados pra um valor de corte em MeV
load('workspace\energia_MeV')

FL0C0M0 = [];
FL0C1M0 = [];
FL0C2M0 = [];
FL0C3M0 = [];

for i=1:50385
    if Ma1(i,1)>150
        FL0C0M0 = [FL0C0M0; L0C0M0(i,:)];
    end
end

for i=1:50385
    if Ma1(i,2)>150
        FL0C1M0 = [FL0C1M0; L0C1M0(i,:)];
    end
end

for i=1:50385
    if Ma1(i,3)>150
        FL0C2M0 = [FL0C2M0; L0C2M0(i,:)];
    end
end

for i=1:50385
    if Ma1(i,4)>150
        FL0C3M0 = [FL0C3M0; L0C3M0(i,:)];
    end
end


tempo = toc/60

stop = 1;

%% Normalizando os dados

NFL0C0M0 = FL0C0M0(:,:);
for i=1:23344
    div = norm(NFL0C0M0(i,:));
for j=1:7
    NFL0C0M0(i,j)=NFL0C0M0(i,j)/div;
end
end

NFL0C1M0 = FL0C1M0(:,:);
for i=1:23634
    div = norm(NFL0C1M0(i,:));
for j=1:7
    NFL0C1M0(i,j)=NFL0C1M0(i,j)/div;
end
end

NFL0C2M0 = FL0C2M0(:,:);
for i=1:42738
    div = norm(NFL0C2M0(i,:));
for j=1:7
    NFL0C2M0(i,j)=NFL0C2M0(i,j)/div;
end
end

NFL0C3M0 = FL0C3M0(:,:);
for i=1:43371
    div = norm(NFL0C3M0(i,:));
for j=1:7
    NFL0C3M0(i,j)=NFL0C3M0(i,j)/div;
end
end

tempo = toc/60

stop = 1;


