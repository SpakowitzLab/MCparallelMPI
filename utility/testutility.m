clear; close all
r = load('r101');

% parameters
boxl = 20;
DEL = 1;
V = 0.1;
lksample = 10;
NP = 2000;
EPS = 0.01; G = 5; N = 8;
L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);
LP=L0/(2*EPS);

% structure factor
[k,S]=scalc(r,boxl,lksample);
figure;loglog(k,S)

% density
[PHIA,PHIB] = r_to_phi(r,boxl,DEL,V);
figure;hist(PHIA+PHIB,50);

% end-to-end distribution
[X,P]=pcalc(r,NP,L0);
figure;hold;
plot(X,P);pplot(N*G,EPS*N*G)