clc, clear all, close all, warning off
addpath(genpath(pwd));

n     = 5000;           % signal dimension 
m     = ceil(.2*n);    % number of measurements
s     = ceil(.002*n);   % signal sparsity
r     = 0.05;          % sign flipping ratio
v     = 0.2;           % correlation parameter
S     = 10;           % number of trials
k=10;
sigma=0.1;

test  = 's_*';           % test \in {'s_*','m/n','r','v','n'}
type  = 'Cor';         % type \in {'Ind', 'Cor'}
switch test
    case 's_*', test0 = ceil(2:4:10);
    case 'm/n', test0 = 0.1:0.2:1.5;
    case 'r',   test0 = 0.05:0.05:0.15;
    case 'v',   test0 = 0.1:0.3:0.7; type = 'Cor'; 
end

algo  = 9;
SNR   = zeros(algo,S,length(test0));
TIME  = zeros(algo,S,length(test0));
HD    = zeros(algo,S,length(test0));
HE    = zeros(algo,S,length(test0));
SPA   = zeros(algo,S,length(test0));

for j = 1:length(test0)
    switch test
        case 's_*',   s=test0(j);
        case 'm/n',   m=ceil(test0(j)*n);
        case 'r',     r=test0(j);
        case 'v',     v=test0(j);    
        case 'n',     n=test0(j); s=ceil(0.00*n); m=ceil(n/2);    
    end
    rand('seed',0);    
    randn('seed',0);  
    snr   = zeros(algo,S);
    time  = zeros(algo,S);
    hd    = zeros(algo,S);
    he    = zeros(algo,S);
    spa   = zeros(algo,S);
    K=s;
    for ii = 1 : S    
        % --------- Generate data ---------------------
%        [A0,c,ctrue,xtrue]...
%              =  random1bcs(type,m,n,s,0.1,r,v);
    seednum = k*ii+1000;
      rand('seed',seednum);   % fix seed 
     randn('seed', seednum); % fix seed 
%    
    xtrue = zeros(n,1);
    rp = randperm(n);
    xtrue(rp(1:K)) = sign(randn(K,1));
    xtrue = xtrue/norm(xtrue);
    supp = find(xtrue);
    SIGMA = zeros(n,n);
    for i = 1:n
        for j = 1:n
           SIGMA(i,j) =v^(abs(i-j));
        end
    end
    Mu = zeros(1,n);
    A0=  mvnrnd(Mu,SIGMA,m);
    ctrue = A0*xtrue;
    noise = sigma*randn(m,1);
    c = sign( ctrue + noise);
    nflip = floor(r*m);
    indxflip = randperm(m);
    indxflip = indxflip(1:nflip);
    c(indxflip) = -c(indxflip);   


         
        x0   = zeros(n,1);
        X    = zeros(n,algo);

        fprintf('------ BIHT  ----------\n')
        in       = 1;
        t        = tic;
        xb = BIHT_1(c, A0, s);
        time(in) = time(in) + toc(t); 
        X(:,in)  = xb/norm(xb); 
        L0       = ceil(nnz(sign(A0*xb)-ctrue));
     
     

        fprintf('------ AOP  ------\n')
        in       = in+1;
        L        = L0; 
        t        = tic;
        xa       = BIHT_AOP_flip(c, A0, x0, s, L, 1, 100, 1);
        time(in) = time(in) + toc(t); 
        X(:,in)  = xa/norm(xa);
        
         fprintf('------ LinProj ----------\n')
         in=in+1;
          t=tic,
         tx = A0'*c/m;
         xl = LinProj(tx,sqrt(s),1,.5,100,50,1e-6);
          time(in) = time(in) + toc(t);
          X(:,in) =xl/norm(xl);    
     
        fprintf('------ PinBalAOP  ------\n')
        in       = in+1;
        t        = tic;
        xpi      = PIHT_AOP_flip(c, A0, x0, s, L, 1, 100, 1, .05); 
        time(in) = time(in) + toc(t); 
        X(:,in)  = xpi/norm(xpi);

        fprintf('------  PDASC  ------\n')
        in       = in+1;
        t        = tic;
        xpd      = pdasc(A0,A0',c);
        time(in) = time(in) + toc(t ); 
        X(:,in)  = xpd/(norm(xpd)+1e-10); 
        fprintf('------ GPSP----------\n') 
        in       = in+1;
        t        = tic;
         k        = ceil(0.01*m); 
        out      = GPSP(A0,c,s,k);        
        time(in) = time(in) + toc(t); 
        X(:,in)  = out.x/norm(out.x);
        fprintf('------ WPDASC ------\n') 
        in       = in+1;
        t        = tic;
        xw       = wpdasc(A0,A0',c);
        time(in) = time(in) + toc(t ); 
        X(:,in)  = xw/(norm(xw)+1e-10); 
       
           fprintf('------ MCP------\n') 
        in       = in+1;
        t        = tic;
        xm=mcp_1bit(c,A0,x0,s);
        time(in) = time(in) + toc(t ); 
        X(:,in)  = xm/(norm(xm)); 
        
       fprintf('------ MCPWP ------\n') 
        in       = in+1;
        t        = tic;
        xmc       = mcp_pdasc(A0,A0',c);
        time(in) = time(in) + toc(t ); 
        X(:,in)  = xmc/(norm(xmc)+1e-10); 
        for alg=1:algo
            snr(alg,ii)= - 20*log10(norm(xtrue-X(:,alg)));    
            hd(alg,ii)=  + nnz(sign(A0*X(:,alg))-c)/m;
            he(alg,ii)=  + nnz(sign(A0*X(:,alg))-ctrue)/m;
            spa(alg,ii)=  + nnz(X(:,alg)==0&xtrue~=0)/s;
        end
    end
    

    SNR(:,:,j) =  snr;
    TIME(:,:,j) =  time;
    HD(:,:,j) = hd;
    HE(:,:,j) = he;
    SPA(:,:,j) = spa;
end
 
 
% algs = {'BIHT','AOP','LinProj','PBAOP','PDASC','GPSP','WPDASC','MCP','MCPDASC',};
algs = {'A1','A2','A3','A4','A5','A6','A7','A8','A9'};
ylab = {'SNR','HD','HE','SMR'};
figure('Renderer', 'painters', 'Position', [800 200 600 400])  
I = [-0.07 -0.02 0.03];
J = [-0.02 -0.02 -0.04];
for item=1:3 
     P=SNR*(item==1)+HD*(item==2)+HE*(item==3)+SPA*(item==4);  
    for i = 1:length(test0)
        sub  = subplot(3,length(test0),(item-1)*3+i);
        pos1 = get(sub, 'Position');  
        if item < 3
            boxplot(P(:,:,i)', 'whisker', 4); set(gca,'xtick',[]); hold on
        else
            boxplot(P(:,:,i)',algs, 'whisker', 4); hold on
        end
        
        if item>1;  set(gca,'ytick',[0,0.1,0.2,0.3,0.4]); end

        if i==1; ylabel(ylab{item},'FontSize',8); end
        if item==1 
           set(0,'DefaultAxesTitleFontWeight','normal');
           title(strcat(test,'=',num2str(test0(i))),'FontSize',8);  
        end
        set(sub, 'Position',pos1+[I(i),J(item),0.06,0.05] )
    end

end


if test == 's_*'; test = 's'; end
if test == 'm/n'; test = 'm'; end
path = strcat('outputs\boxplot-',test);  
saveas(figure(1), strcat(path,'.fig'));  