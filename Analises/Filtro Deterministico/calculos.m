load('dados\Lado-Canal.mat');
load('workspace\energia_MeV.mat');

%% Criando matrizes separando por lado, canal e modulo
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


