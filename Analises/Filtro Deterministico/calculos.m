load('dados\Lado-Canal.mat');
load('workspace\energia_MeV.mat');

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


