clear;close all
% inputs:
% r101, an example file for bead coordinates
% pdata, analytical solution for end-to-end distribution
%   pdata available shifan@tower12:~/testdir/pdata

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
Ree = 2;

%% Single chain conformation tests
% end-to-end distribution
[X,P]=pcalc(r,NP,L0);
figure;hold;
plot(X,P);pplot(N*G,EPS*N*G)
xlabel('R/L');ylabel('P(R/L)')

% bead-bead distribution
R = rcalc(r,NP);
NK = (0:N*G-1)'/G;
NM = EPS*G;
DM = EPS*NK*G;
RM = NM-(1/2)*(1-exp(-2*NM));
RD = DM-(1/2)*(1-exp(-2*DM));

figure;hold
plot(NK,sqrt(R)/Ree,'o',NK,sqrt(RD/RM),'-')
set(gca,'xscale','log');set(gca,'yscale','log')
xlabel('J/G');ylabel('R_J/R_M')

%% Collective conformation tests
% density
[PHIA,PHIB] = r_to_phi(r,boxl,DEL,V);
figure;hist(PHIA+PHIB,50);
xlabel('\phi_A+phi_B');ylabel('P(\phi)')

% radial distrbution function
[R,g] = gcalc(r,boxl);
figure;
plot(R,g(:,1),R,g(:,2),R,g(:,3),R,g(:,4));
xlabel('r');ylabel('g(r)')

% structure factor
[k,S]=scalc(r,boxl,lksample);
figure;loglog(k,S)
xlabel('q');ylabel('S(q)')
