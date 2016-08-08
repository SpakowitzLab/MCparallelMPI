clear;close all

% OPTIONS
PLOTSIM = 1;    % plot simulation structure factor

% plot parameters
NREP = [1:79];  % number of replicas
NSNAP = 50;  % snapshots to average

% simulation parameters
EPS = 0.01;  % inter-bead segment rigidity (in unit of 2lp)
LAM = 0.0;     % degree of chemical correlation

% simulation constants
boxl = 20;   % edge size of simulation
Ree = 2.0;   % average end-to-end distance of a monomer
G = 5;       % number of beads per monomer
N = 8;
FA = 0.5;    % chemical fraction of A species
NP = 2000;   % number of polymers
NK = (0:N*G-1)'/G;
L0=2*EPS*(-.5+.5*exp(-2*EPS*G)+EPS*G)^(-.5);

% add/define paths
%loaddir = '../data/';    % data directory
[pathstr,name,ext] = fileparts(pwd);
title = 'randcopoly';ind = findstr(pathstr, title);
loaddir = strcat('../../../sim-',pathstr(ind:end),'/data/');
savedir = 'savedata/';   % save directory
addpath('misc/');
addpath('../utility/');

% load CHI parameters
CHI = load(strcat(loaddir,'chi'));  % adjusted CHI values
CHI = CHI(end,2:end);

figure(1);hold;set(gca,'fontsize',20);
figure(2);hold;set(gca,'fontsize',20);

if (PLOTSIM)
    % plot simulation results
    for REP=NREP
      if max(NREP) == 1
        col = 1;
      else
        col = (REP-1)/(max(NREP)-1);
      end
      Savg = [];pdfavg = [];
      for SNAP=NSNAP
          fprintf('REP = %d, SNAP = %d\n',REP,SNAP)
          r=dlmread(strcat(loaddir,sprintf('r%dv%d',SNAP,REP)));

          % calculate bead-to-bead distance
          [cosTheta,pdf,n,S]=tcalc(r,NP,L0);
          Savg = [Savg,S];
          pdfavg = [pdfavg,pdf'];
      end
      pdfavg = mean(pdfavg,2);
      Savg = mean(Savg);

      figure(1);plot(cosTheta,pdfavg,'o-','linewidth',1.5,'color',[col 0 1-col]);
      figure(2);plot(CHI(REP)*G,Savg,'o','linewidth',1.5,'color',[col 0 1-col]);
    end
end

figure(1);
xlabel('cos(\theta)');ylabel('P(cos(\theta))');box on

figure(2);
xlabel('\chiG');ylabel('<S>');box on
