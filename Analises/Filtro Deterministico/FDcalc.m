%clear all; close all; clc;

load('muonMa1.mat');
load('noiseMa1.mat');

IdSinal = muonMa1(:,1) + muonMa1(:,2) + muonMa1(:,3) + muonMa1(:,4);
IdRuido = noiseMa1(:,1) + noiseMa1(:,2) + noiseMa1(:,3) + noiseMa1(:,4);

%% ROC Filtro Casado

pmin = -1020;
pmax = 3020;
pontos = 2000;

psoma = (pmax+abs(pmin))/pontos;
patamar = pmin;
PD = zeros(pontos);
FA = zeros(pontos);
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
    PD(i) = pd; % preenchendo o vetor
    FA(i) = fa;
    pd = 0;
    fa = 0;
    patamar = patamar + psoma;
end

hold on
plot(FA, PD, '-xr')
legend('Filtro Estocastico (simplificado)', 'Filtro Estocastico (completo)', 'Filtro Casado');
grid
% title('ROC')
% xlabel('% FA')
% ylabel('% PD')