function R = rcalc(r,NP)

NB=length(r);
N=NB/NP;

%%%% start calculations %%%%
R = zeros(N,1);
for IK=1:NP
  K1=(IK-1)*N+1;
  K2=[1:N]+(IK-1)*N;
  R=R+sum(power(r(K2,1:3)-repmat(r(K1,1:3),N,1),2),2);
end
R = R/NP;

end
