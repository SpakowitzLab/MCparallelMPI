clear;
%close all
addpath('misc/');
addpath('../utility/');

% parse replicas
%dir = '../data/';  % data directory
%dir = '../../../quinn/MCPoly/MCparallelMPI/lamNeg1_r7_27_16/';
%cof = load(strcat(dir,'cofData'));
chi = load(strcat(dir,'chi'));
node = load(strcat(dir,'nodeNumber'));
[REPS,CHI,NTIME,NREP] = replica(chi,node);

% simulation parameters
G = 5;

figure;hold
for ii = 1:NREP
 plot(1:10:NTIME,CHI(1:10:NTIME,ii)*G,'linewidth',1.5)
end

figure;hold
for ii = 250:1:290
    plot(CHI(250,:)*G,CHI(ii,:)*G,'rx')
end
