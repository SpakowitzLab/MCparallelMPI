function [R,g] = gcalc(r,boxl)
% INPUTS:
% OUTPUTS: 
% R, radial separation
% G, radial distribution
% G = [GAA,GAB,GBB,GTOT]

%% Pre-process data
% parameters
NB = length(r);     % total number of beads
RHO = NB/(boxl^3);  % bead density

% NTOT = NB*100; % total # of sampling pairs
NTOT = 2e6;
RMAX = boxl/4;      % sampling range
DEL = 0.02;         % sampling step size
N = round(RMAX/DEL);% number of sampling steps

% apply periodic boundary conditions to coordinates
RB = r;
for K=1:NB
   RB(K,1)=RB(K,1)-1.*round(RB(K,1)/boxl-0.5)*boxl;
   RB(K,2)=RB(K,2)-1.*round(RB(K,2)/boxl-0.5)*boxl;
   RB(K,3)=RB(K,3)-1.*round(RB(K,3)/boxl-0.5)*boxl;
end

R = zeros(N,1);
g = zeros(N,4);
for J = 1:N
    R(J) = DEL/2+(J-1)*DEL;
end

for I = 1:NTOT
   I1=ceil(rand*NB);
   I2=ceil(rand*NB);

   R1(1:3)=RB(I1,1:3);
   DR(1)=boxl*round((RB(I2,1)-RB(I1,1))/boxl);
   DR(2)=boxl*round((RB(I2,2)-RB(I1,2))/boxl);
   DR(3)=boxl*round((RB(I2,3)-RB(I1,3))/boxl);
   R2(1:3)=RB(I2,1:3)-DR(1:3);

   DIST = sqrt((R1(1)-R2(1))^2+(R1(2)-R2(2))^2+(R1(3)-R2(3))^2);
   IND=round(3-RB(I1,4)-RB(I2,4));
   RIND=round((DIST+DEL/2)/DEL)+1;

   if ((RIND<N) && (I1~=I2))
      g(RIND,IND)=g(RIND,IND)+1;
      % Use Fredrickson PRL 1991 ORDER PARAMETER
      g(RIND,4)=g(RIND,4)+1;
   end
end

% %normalization from JJDepablo
for J = 1:N
   g(J,1) = NB*g(J,1)/(4*pi*RHO*NTOT*DEL);
   g(J,2) = NB*g(J,2)/(4*pi*RHO*NTOT*DEL);
   g(J,2) = g(J,2)/2;
   g(J,3) = NB*g(J,3)/(4*pi*RHO*NTOT*DEL);
   g(J,4) = NB*g(J,4)/(4*pi*RHO*NTOT*DEL);
   
   g(J,1)=g(J,1)/(R(J)*R(J));
   g(J,2)=g(J,2)/(R(J)*R(J));
   g(J,3)=g(J,3)/(R(J)*R(J));
   g(J,4)=g(J,4)/(R(J)*R(J));
end

end
