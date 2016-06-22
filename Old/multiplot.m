noIntBind=zeros(8,2);
intBind=zeros(8,2);
for j=1:8
    figure(j)
    
    disp('reading params')
    params=dlmread(['bigData/con' num2str(j) '/data/params']); % parameters
    NT=params(1);
    N=params(2);
    NB=params(3);
    NP=params(4);
    disp('NP')
    disp(NP)
    NT=params(5);
    G=params(6);

    L0=params(7);
    CHI0=params(8);
    Fpoly=params(9);
    LBox=params(10);
    EU = params(11)
    EM = params(12)
    HP1_Bind = params(13)
    LK = params(14);
    FRMMETH = params(15);
    F_METH = params(16);
    LAM_METH = params(17);

    subplot(1,2,1)
    disp('csv to matlab...')
    data=dlmread(['bigData/con' num2str(j) '/data/r499']);

    [rows, cols]=size(data);

    disp('Applying periodic condition')
    data(:,1:3)=mod(data(:,1:3),LBox);


    disp('coloring...')
    c=zeros(rows,3);
    s=ones(rows,1)*30;
    for row=1:rows
        if data(row,4) % one is read blue is zero
            c(row,:)=[1,0,0]; %red
        else
            c(row,:)=[0,0,1]; %blue
        end
    end

    disp('plotting ...')
    scatter3(data(:,1),data(:,2),data(:,3),s,c,'filled')
    title(['without interaction mu=' num2str(-10+2*j)])
    
    
    subplot(1,2,2)
    disp('csv to matlab...')
    data=dlmread(['bigData/con' num2str(j) '/data/r999']);

    [rows, cols]=size(data);

    disp('Applying periodic condition')
    data(:,1:3)=mod(data(:,1:3),LBox);


    disp('coloring...')
    c=zeros(rows,3);
    s=ones(rows,1)*30;
    for row=1:rows
        if data(row,4) % one is read blue is zero
            c(row,:)=[1,0,0]; %red
        else
            c(row,:)=[0,0,1]; %blue
        end
    end

    disp('plotting ...')
    scatter3(data(:,1),data(:,2),data(:,3),s,c,'filled')  
    title('with interaction')
    
    data=dlmread(['bigData/con' num2str(j) '/data/r499']);
    u0=dlmread(['bigData/con' num2str(j) '/data/u0']); % u0
    METH=logical(u0(:,5));
    nmeth=sum(METH);
    numeth=length(METH)-nmeth;
    
    msum=sum(METH.*data(:,4));
    wsum=sum((~METH).*data(:,4));
    noIntBind(j,1)=msum/nmeth;
    noIntBind(j,2)=wsum/numeth;

    data=dlmread(['bigData/con' num2str(j) '/data/r999']);
    msum=sum(METH.*data(:,4));
    wsum=sum((~METH).*data(:,4));
    intBind(j,1)=msum/nmeth;
    intBind(j,2)=wsum/numeth;    

    xlim([0 LBox])
    zlim([0 LBox])
    ylim([0 LBox])
    
end
figure(9)
hold on
plot(-8:2:6,noIntBind(:,1),'-or')
plot(-8:2:6,noIntBind(:,2),'-ob')
plot(-8:2:6,intBind(:,1),'-xr')
plot(-8:2:6,intBind(:,2),'-xb')
legend('no interaction, methalated',...
       'no interaction, unmethalated',...
       'interaction, unmethalated',...
       'interaction, methalated')
ylabel('fraction bound')
xlabel('x')



