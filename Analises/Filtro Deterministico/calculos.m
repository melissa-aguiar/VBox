clear all
close all
clc


%% normalizando os dados que foram filtrados pra >150MeV
load('dados filtrados Lado-Canal/FL0C0.mat');
load('dados filtrados Lado-Canal/FL0C1.mat');
load('dados filtrados Lado-Canal/FL0C2.mat');
load('dados filtrados Lado-Canal/FL0C3.mat');
load('dados filtrados Lado-Canal/FL1C0.mat');
load('dados filtrados Lado-Canal/FL1C1.mat');
load('dados filtrados Lado-Canal/FL1C2.mat');
load('dados filtrados Lado-Canal/FL1C3.mat');

NFL0C0 = FL0C0(:,(4:10));
NFL0C1 = FL0C1(:,(4:10));
NFL0C2 = FL0C2(:,(4:10));
NFL0C3 = FL0C3(:,(4:10));
NFL1C0 = FL1C0(:,(4:10));
NFL1C1 = FL1C1(:,(4:10));
NFL1C2 = FL1C2(:,(4:10));
NFL1C3 = FL1C3(:,(4:10));

for i=1:1480216
    div = norm(NFL0C0(i,:));
for j=1:7
    NFL0C0(i,j)=NFL0C0(i,j)/div;
end
end

for i=1:1478495
    div = norm(NFL0C1(i,:));
for j=1:7
    NFL0C1(i,j)=NFL0C1(i,j)/div;
end
end

for i=1:2709185
    div = norm(NFL0C2(i,:));
for j=1:7
    NFL0C2(i,j)=NFL0C2(i,j)/div;
end
end

for i=1:2685094
    div = norm(NFL0C3(i,:));
for j=1:7
    NFL0C3(i,j)=NFL0C3(i,j)/div;
end
end

for i=1:1508572
    div = norm(NFL1C0(i,:));
for j=1:7
    NFL1C0(i,j)=NFL1C0(i,j)/div;
end
end

for i=1:1539631
    div = norm(NFL1C1(i,:));
for j=1:7
    NFL1C1(i,j)=NFL1C1(i,j)/div;
end
end

for i=1:2651738
    div = norm(NFL1C2(i,:));
for j=1:7
    NFL1C2(i,j)=NFL1C2(i,j)/div;
end
end

for i=1:2693031
    div = norm(NFL1C3(i,:));
for j=1:7
    NFL1C3(i,j)=NFL1C3(i,j)/div;
end
end

stop = 1


%% Filtrando os dados pra um valor de corte em MeV
tic

G1L0C0 = [];
aux1 = [];

j=0;

for i=1:3203840
    j=j+1;
    if Ma(i,1)>1000 %muda o canal muda o Ma
        aux1 = [aux1; L0C2(i,:)];
    end
    
    if j==40000
        G1L0C0 = [G1L0C0; aux1];
        aux1 = [];
        j=0;
    end
        
  end
G1L0C0 = [G1L0C0; aux1];

tempo = toc/60