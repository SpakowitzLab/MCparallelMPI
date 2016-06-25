% plot structure factors
clear; close all
addpath('../utility/');

% simulation parameters
boxl = 20;
DEL = 1;
V = 0.1;
Ree = 2;
EPS = 0.01;
LAM = -0.75;
FA = 0.5;
N = 8;
G = 5;
KAP = 10;
CHI = load('../data/cof');

% plot parameters
lksample = 20;
NREP = 1:4:40;  % number of replicas
NAVG = 10:15;  % snapshots to average

% results to return
MU = [];
for REP=NREP
    REP
    mu = [];
    for AVG=NAVG
        r=dlmread(sprintf('../data/r%dv%d',AVG,REP));
        mu=[mu,widom(r,LAM,FA,N,G,CHI(REP),KAP,boxl,DEL,V)];
    end
    MU=[MU,mean(mu)];
end

% plot simulation results
figure;plot(1:length(MU),MU,'o-')