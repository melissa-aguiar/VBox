% clear all
% close all
% clc

% load('NFL0C0M0');
% load('NFL0C1M0');
% load('NFL0C2M0');
% load('NFL0C3M0');

%% Plotando os pulsos Lado0-Modulo0-Canais0123

figure
plot(1:7,NFL0C0M13(:,:))
title('Lado 0 Canal 0 Modulo 13')
grid on

figure
plot(1:7,NFL0C1M13(:,:))
title('Lado 0 Canal 1 Modulo 13')
grid on

figure
plot(1:7,NFL0C2M13(:,:))
title('Lado 0 Canal 2 Modulo 13')
grid on

figure
plot(1:7,NFL0C3M13(:,:))
title('Lado 0 Canal 3 Modulo 13')
grid on

%% Calculando e plotando autovalores e autovetores

[COEFF0, SCORE0, LATENT0] = pca(FL0C0M13);
[COEFF1, SCORE1, LATENT1] = pca(FL0C1M13);
[COEFF2, SCORE2, LATENT2] = pca(FL0C2M13);
[COEFF3, SCORE3, LATENT3] = pca(FL0C3M13);

figure
subplot(2,2,1);
hist(COEFF0)
title('COEFF Lado 0 Canal 0 Modulo 13')
grid
subplot(2,2,2);
hist(COEFF1)
title('COEFF Lado 0 Canal 1 Modulo 13')
grid
subplot(2,2,3);
hist(COEFF2)
title('COEFF Lado 0 Canal 2 Modulo 13')
grid
subplot(2,2,4);
hist(COEFF3)
title('COEFF Lado 0 Canal 3 Modulo 13')
grid

figure
subplot(2,2,1);
plot(LATENT0,'-x')
title('LATENT Lado 0 Canal 0 Modulo 0')
grid
subplot(2,2,2);
plot(LATENT1,'-x')
title('LATENT Lado 0 Canal 1 Modulo 0')
grid
subplot(2,2,3);
plot(LATENT2,'-x')
title('LATENT Lado 0 Canal 2 Modulo 0')
grid
subplot(2,2,4);
plot(LATENT3,'-x')
title('LATENT Lado 0 Canal 3 Modulo 0')
grid
