% plot structure factors
clear; close all
addpath('../utility/');

% simulation parameters
boxl = 20;
Ree = 2;
EPS = 0.01;
LAM = -0.75;
G = 5;

% plot parameters
lksample = 20;
NREP = [1:2:10,19,29,39];  % number of replicas
NAVG = 10:1:15;  % snapshots to average

figure;hold
% plot simulation results
for REP=NREP
    col = (REP-1)/(max(NREP)-1)
    savg = [];
    for AVG=NAVG
        r=dlmread(sprintf('../data/r%dv%d',AVG,REP));
        [k,s]=scalc(r,boxl,lksample);
        savg = [savg,s];
    end
    plot(k*Ree,mean(savg,2),'color',[col 0 1-col]);
end

% plot mean-field theory
CHI = load('../data/cof');
filename = sprintf('../utility/sdata/Seps%.3flam%.2f',EPS,LAM);
spinodal = load('../utility/sdata/chivals');
CHIS = spinodal(spinodal(:,1)==EPS & spinodal(:,2)==LAM,3);
S = load(filename);
KS = S(:,1);
for REP=NREP
    col = (REP-1)/(max(NREP)-1);
    
    if CHI(REP) < CHIS/G
        SS = 1./(-2*CHI(REP)+1./S(:,2));
        plot(KS,SS,'--','color',[col 0 1-col]);
    end
end

% plot processing
set(gca,'xscale','log');set(gca,'yscale','log')