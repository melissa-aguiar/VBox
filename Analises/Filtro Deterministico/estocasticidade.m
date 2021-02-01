clear all
close all
clc

tic

load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL0C0.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL0C1.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL0C2.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL0C3.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL1C0.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL1C1.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL1C2.mat');
load('C:\Users\melis\Desktop\VBox\Analises\Filtro Deterministico\dados normalizados e filtrados Lado-Canal\NFL1C3.mat');



%SmpMuon = load('Sample_Muons.txt'); % | LADO[2] | MODULO[64] | CANAL[4] | AMOSTRAS[7] | ... CANAL 0, 1 é uma CELULA; 2 e 3 é outra

%% LADO 0 - TODOS CANAIS
figure
plot(1:7,NFL0C0(1:100000,:))
title('Lado 0 Canal 0')
grid on

tempo = toc/60

stop = 1

figure
plot(1:7,NFL0C1(1:100000,:))
title('Lado 0 Canal 1')
grid on

figure
plot(1:7,NFL0C2(1:100000,:))
title('Lado 0 Canal 2')
grid on

figure
plot(1:7,NFL0C3(1:100000,:))
title('Lado 0 Canal 3')
grid on



tempo = toc/60

stop = 1

figure
for i=1:1478495
    plot(NFL0C1(i,:))
    hold on
end
title('Lado 0 Canal 1')
grid

figure
for i=1:2709185
    plot(NFL0C2(i,:))
    hold on
end
title('Lado 0 Canal 2')
grid

figure
for i=1:2685094
    plot(NFL0C3(i,:))
    hold on
end
title('Lado 0 Canal 3')
grid



%% PCA
[COEFF0, SCORE0, LATENT0] = pca(c0n);
[COEFF1, SCORE1, LATENT1] = pca(c1n);
[COEFF2, SCORE2, LATENT2] = pca(c2n);
[COEFF3, SCORE3, LATENT3] = pca(c3n);

figure
subplot(2,2,1);
hist(COEFF0)
title('COEFF SmpMuon Canal 0')
grid
subplot(2,2,2);
hist(COEFF1)
title('COEFF SmpMuon Canal 1')
grid
subplot(2,2,3);
hist(COEFF2)
title('COEFF SmpMuon Canal 2')
grid
subplot(2,2,4);
hist(COEFF3)
title('COEFF SmpMuon Canal 3')
grid

figure
subplot(2,2,1);
plot(LATENT0,'-x')
title('LATENT SmpMuon Canal 0')
grid
subplot(2,2,2);
plot(LATENT1,'-x')
title('LATENT SmpMuon Canal 1')
grid
subplot(2,2,3);
plot(LATENT2,'-x')
title('LATENT SmpMuon Canal 2')
grid
subplot(2,2,4);
plot(LATENT3,'-x')
title('LATENT SmpMuon Canal 3')
grid
