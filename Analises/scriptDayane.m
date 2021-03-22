clear all
close all
clc

SmpMuon = load('Sample_Muons.txt'); % | LADO | MODULO | CELULA | AMOSTRAS |

ped = 0;

%% gera os dados
ruido = ped + randn(20000,7); % um pedestal + distribuição gaussiana
ruidoDes = ruido(1:10000,:); %depois usa no branqueamento
ruidoTes = ruido(10001:end,:); %depois usa na aplicação do filtro estocastico pra calcular Id e Ir

ampTrue = zeros(1,10000);
sinalTes = zeros(10000,7);

for i=1:10000 % faz a composição do sinal 10k
    ampTrue(i) = 1023*rand(1); %uma amplitude pra cada evento
    sinalTes(i,:) = ped + randn(1,7) + ampTrue(i)*pegaPulseJitter(); %ruido + distribuiçao uniforme + estocasticidade
end

%% script para branqueamento -- descorrelacionando o ruido
c = cov(ruidoDes); %ruidoDes é descorrelacionado
[V,D] = eig(c); % autovalores e autovetores
W = D^(-.5)*V'; % branqueamento
%ruido branqueado

%% script para PCAs
sinalPuro = zeros(10000,7);
for i=1:10000
    sinalPuro(i,:) = pegaPulseJitter();
end
[COEFF, SCORE, LATENT] = pca((sinalPuro*W')); % matlab nao reconheceu princomp (principal component analysis)
%essa funçao adquire os autovetores e autovalores do sinal

%% alguns parametros para o filtro estocastico
variancia = var(ruidoDes(:,4)); %parametro usado no filtro estocastico
sm = [0.0000 0.0172 0.4524 1.0000 0.5633 0.1493 0.0424]; %sinal medio pra calcular a parte deterministica
mFC = sm*W'; % sinal medio para calcular Id

%% aplicacao do filtro estocastico
for pc=7:7 % 7 é randomico? é por causa do numero de amostras?

    N = pc;  % numero de PCAs
    mEstimacao = mFC*COEFF(:,1:N);
    
    rRuido = ((ruidoTes-ped)*W')*COEFF(:,1:N); %sinal de entrada do ruido
    rSinal = ((sinalTes-ped)*W')*COEFF(:,1:N);
   
    No = variancia*2;
    lambda = LATENT;
    h1 = zeros(7,7); % vai ser a parte constante na formula de IR
    h2 = zeros(7,7); % vai ser a parte constante na formula de ID
    for i=1:N % de 1 ate o numero de pca
        h1 = h1 + ((lambda(i))./((lambda(i))+variancia))*(COEFF(:,i)*COEFF(:,i)');
        h2 = h2 + ((1./((lambda(i))+variancia)))*(COEFF(:,i)*COEFF(:,i)');
    end
    
    %------------------ acha a parte estocastica --------------------------
    IrRuido = zeros(size(ruidoTes,1),1); 
    IrSinal = zeros(size(sinalTes,1),1);
    for ev=1:size(ruidoTes,1) % testa qual o valor do evento usando a formula inteira do artigo
        IrRuido(ev) = (1/No)*((rRuido(ev,:)*COEFF(:,1:N)')*h1*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IrSinal(ev) = (1/No)*((rSinal(ev,:)*COEFF(:,1:N)')*h1*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

    %------------------ acha a parte deterministica -----------------------
    IdSinal = zeros(size(sinalTes,1),1);
    IdRuido = zeros(size(ruidoTes,1),1);
    for ev=1:size(ruidoTes,1)
        IdRuido(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IdSinal(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

    
    FCestRuido = IdRuido + IrRuido; % saida do filtro pra deteccao para o ruido
    FCestSinal = IdSinal + IrSinal; % saida do filtro pra deteccao para o sinal

    %% estimacao da amplitude 
    b1 = COEFF(:,1:N)*COEFF(:,1:N)';
    b2 = (1/No)*(COEFF(:,1:N)*h1*COEFF(:,1:N)');
    b3 = (mEstimacao*COEFF(:,1:N)')*h2*COEFF(:,1:N)';
    
    ampRuido = zeros(size(ruidoTes,1),1);
    ampSinal = zeros(size(sinalTes,1),1);
    a = (1/No)*((mEstimacao*COEFF(:,1:N)')*h1*(mEstimacao*COEFF(:,1:N)')');
    b = (mEstimacao*COEFF(:,1:N)')*h2*(mEstimacao*COEFF(:,1:N)')';
    cs=0;
    cr=0;
    for i=1:size(sinalTes,1)
        ra = b*b+4*a*FCestSinal(i);
        if ra<0
            ra=0;
            cs=cs+1;
        end
        ampSinal(i) = (-b+sqrt(ra))/(2*a); % amplitude do sinal usando a saida do filtro casado
    end
    for i=1:size(ruidoTes,1)
        ra = b*b+4*a*FCestRuido(i);
        if ra<0
            ra=0;
            cr=cr+1;
        end
        ampRuido(i) = (-b+sqrt(ra))/(2*a); % amplitude do ruido usando a saida do filtro casado
    end
   
end


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

load('Filtro Deterministico\muonMa1.mat');
load('Filtro Deterministico\noiseMa1.mat');

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

% hold on
% plot(FA1, PD1, '-mx');
% 
% % legend('Filtro Casado', 'Filtro Estocástico (Célula)','Filtro Estocástico (Canal)','Filtro Estocástico (Módulo)')
% title('ROC - Results')
% legend('Matched Filter', 'Stochastic Filter (Cell)', 'Stochastic Filter (Channel)', 'Stochastic Filter (Module)');
% xlabel('% FA')
% ylabel('% PD')
% grid
% stop = 1

%%
figure

plot(FA1, PD1, '-b*');

hold on
plot(FA2, PD2, '-rs');

hold on
plot(FA3, PD3, '-.g', 'LineWidth', 2);

hold on
plot(FA, PD, '-mx');


grid
title('ROC - Module')
legend('Stochastic Filter', 'Deterministic Component', 'Matched Filter', 'Stochastic Component');
xlabel('% FA')
ylabel('% PD')
