function [PHIA,PHIB] = r_to_phi(r,LBOX,DEL,V)
%% calculate PM0 densities from bead positions
% LBOX = 20; DEL = 1; V = 0.1;

%%
R = r(:,1:3);
AB = r(:,4);
NB = length(r);

NBINX=round(LBOX/DEL);
NBIN = NBINX^3;

PHIA = zeros(NBIN,1);
PHIB = zeros(NBIN,1);
RBIN = zeros(3,1);

%     Cycle through the beads
      
IB=1;
for I=1:NB
    RBIN(1)=R(IB,1)-round(R(IB,1)/LBOX-0.5)*LBOX;
    RBIN(2)=R(IB,2)-round(R(IB,2)/LBOX-0.5)*LBOX;
    RBIN(3)=R(IB,3)-round(R(IB,3)/LBOX-0.5)*LBOX;

    IX(1)=round(RBIN(1)/DEL+0.5);
    IY(1)=round(RBIN(2)/DEL+0.5);
    IZ(1)=round(RBIN(3)/DEL+0.5);
            
    IX(2)=IX(1)-1;
    IY(2)=IY(1)-1;
    IZ(2)=IZ(1)-1;

    %     Calculate the bin weighting

    WX(2)=(RBIN(1)-IX(1)*DEL)/(IX(2)*DEL-IX(1)*DEL);
    WX(1)=(IX(2)*DEL-RBIN(1))/(IX(2)*DEL-IX(1)*DEL);
    WY(2)=(RBIN(2)-IY(1)*DEL)/(IY(2)*DEL-IY(1)*DEL);
    WY(1)=(IY(2)*DEL-RBIN(2))/(IY(2)*DEL-IY(1)*DEL);
    WZ(2)=(RBIN(3)-IZ(1)*DEL)/(IZ(2)*DEL-IZ(1)*DEL);
    WZ(1)=(IZ(2)*DEL-RBIN(3))/(IZ(2)*DEL-IZ(1)*DEL);

    IX(1)=IX(1)-floor((IX(1)-1)/NBINX) * NBINX;
    IX(2)=IX(2)-floor((IX(2)-1)/NBINX) * NBINX;
    IY(1)=IY(1)-floor((IY(1)-1)/NBINX) * NBINX;
    IY(2)=IY(2)-floor((IY(2)-1)/NBINX) * NBINX;
    IZ(1)=IZ(1)-floor((IZ(1)-1)/NBINX) * NBINX;
    IZ(2)=IZ(2)-floor((IZ(2)-1)/NBINX) * NBINX;

    %     Add volume fraction with weighting to each bin

    for ISX=1:2
      for ISY=1:2
	for ISZ=1:2
          WTOT=WX(ISX)*WY(ISY)*WZ(ISZ);
          INDBIN=IX(ISX)+(IY(ISY)-1)*NBINX+(IZ(ISZ)-1)*NBINX^2;
          if (AB(IB)==1)
    	    PHIA(INDBIN)=PHIA(INDBIN)+WTOT*V/DEL^3;
          else
            PHIB(INDBIN)=PHIB(INDBIN)+WTOT*V/DEL^3;
          end
	end
       end
    end
    IB=IB+1;
end