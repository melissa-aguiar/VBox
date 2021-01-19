clear all
close all
clc

SmpMuon = load('Sample_Muons.txt');

ped = 0;

%% gera os dados
ruido = ped + randn(20000,7);
ruidoDes = ruido(1:10000,:);
ruidoTes = ruido(10001:end,:);

ampTrue = zeros(1,10000);
sinalTes = zeros(10000,7);
for i=1:10000
    ampTrue(i) = 1023*rand(1);
    sinalTes(i,:) = ped + randn(1,7) + ampTrue(i)*pegaPulseJitter();
end

%% script para branqueamento
c = cov(ruidoDes);
[V,D] = eig(c);
W = D^(-.5)*V'; % branqueamento

%% script para PCAs
sinalPuro = zeros(10000,7);
for i=1:10000
    sinalPuro(i,:) = pegaPulseJitter();
end
[COEFF, SCORE, LATENT] = pca((sinalPuro*W')); %%%%%%%%%%%%%matlab nao reconheceu princomp (principal component analysis)

%% alguns parametros para o filtro estocastico
variancia = var(ruidoDes(:,4)); %parametro usado no filtro estocastico
sm = [0.0000 0.0172 0.4524 1.0000 0.5633 0.1493 0.0424];
mFC = sm*W'; % sinal medio para calcular Id

%% aplicacao do filtro estocastico
for pc=7:7

    N = pc;  % numero de PCAs
    mEstimacao = mFC*COEFF(:,1:N);
    
    rRuido = ((ruidoTes-ped)*W')*COEFF(:,1:N);
    rSinal = ((sinalTes-ped)*W')*COEFF(:,1:N);
   
    No = variancia*2;
    lambda = LATENT;
    h1 = zeros(7,7);
    h2 = zeros(7,7);
    for i=1:N
        h1 = h1 + ((lambda(i))./((lambda(i))+variancia))*(COEFF(:,i)*COEFF(:,i)');
        h2 = h2 + ((1./((lambda(i))+variancia)))*(COEFF(:,i)*COEFF(:,i)');
    end
    
    IrRuido = zeros(size(ruidoTes,1),1);
    IrSinal = zeros(size(sinalTes,1),1);
    for ev=1:size(ruidoTes,1)
        IrRuido(ev) = (1/No)*((rRuido(ev,:)*COEFF(:,1:N)')*h1*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IrSinal(ev) = (1/No)*((rSinal(ev,:)*COEFF(:,1:N)')*h1*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

    IdSinal = zeros(size(sinalTes,1),1);
    IdRuido = zeros(size(ruidoTes,1),1);
    for ev=1:size(ruidoTes,1)
        IdRuido(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rRuido(ev,:)*COEFF(:,1:N)')');
    end

    for ev=1:size(sinalTes,1)
        IdSinal(ev) = ((mEstimacao*COEFF(:,1:N)')*h2*(rSinal(ev,:)*COEFF(:,1:N)')');
    end

    FCestRuido = IdRuido + IrRuido; % saida do filtro pra deteccao para sinal
    FCestSinal = IdSinal + IrSinal; % saida do filtro pra deteccao para ruido

    % estimacao da amplitude 
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

stem(ampSinal(1:50));
figure;
stem(ampRuido(1:50));