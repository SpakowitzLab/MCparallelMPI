function mu=widom(r,LAM,FA,N,G,CHI,KAP,LBOX,DEL,V)
% use widom insertion to determine excess free energy
NI = 100;  % number of insertions
NB = length(r); % total number of beads
NP = NB/N/G;    % number of polymers

% read data
R = r(:,1:3);
% calcuate initial chemical potential
[PHIA0,PHIB0]=r_to_phi(r,LBOX,DEL,V);
ECHI0=(DEL^3.)*CHI/V*sum(PHIA0.*PHIB0);
EKAP0=(DEL^3.)*KAP/V*sum((PHIA0+PHIB0-1).^2);

delU=zeros(NI,1);
for INI=1:NI
    % initialize chemical sequence of ghost mol.
    AB=initchem(LAM,FA,N,G);

    % initialize chain configuration of ghost mol.
    Rg = rand(1,3)*LBOX;
    
    % grow the chain
    for I = 2:N*G
        % all possible configurations
        I1 = (I-1):N*G:NB;
        I2 = I:N*G:NB;
        Rall = R(I2,1:3)-R(I1,1:3);

        pick=randi([1,NP],1);
        %Rgrow = Rg(I-1,1:3)+Rall(pick,1:3);
        vrand=randn(1,3);
        Rgrow = Rg(I-1,1:3)+norm(Rall(pick,1:3))*vrand/norm(vrand);
        Rg = [Rg;Rgrow];
    end
    Rghost = [Rg,AB];
    rnew = [r;Rghost];

    % calculate change in chemical potential
    %[PHIAf,PHIBf]=r_to_phi(datanew,N*G,NP+1);
    [PHIAf,PHIBf]=r_to_phi(rnew,LBOX,DEL,V);
    ECHIf=(DEL^3.)*CHI/V*sum(PHIAf.*PHIBf);
    EKAPf=(DEL^3.)*KAP/V*sum((PHIAf+PHIBf-1).^2);
    
    delU(INI) = (ECHIf-ECHI0)+(EKAPf-EKAP0);
end
mu=-log(mean(exp(-delU)));

end

function AB=initchem(LAM,FA,N,G)
PAA=FA*(1.-LAM)+LAM;
PBB=FA*(LAM-1.)+1.;
PAB=1.-PBB;

AB=zeros(N*G,1);

% start generating polymer chemistry
IB=1;
TEST=rand;
if (TEST<FA)
    AB(IB)=1;
else
    AB(IB)=0;
end

IB=IB+1;

for K=2:G
    AB(IB)=AB(IB-1);
    IB=IB+1;
end

for J=2:N
    TEST=rand;
    if (AB(IB-1)==1)
        if (TEST<=PAA)
            AB(IB)=1;
        else
            AB(IB)=0;
        end
    else
        if (TEST<=PAB)
            AB(IB)=1;
        else
            AB(IB)=0;
        end
    end
    IB=IB+1;

    for K=2:G
        AB(IB)=AB(IB-1);
        IB=IB+1;
    end
end

end