clc
clear all 
close all
warning off
addpath(genpath(fileparts(mfilename('fullpath'))));
% All available solvers:
FIELDS={'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP','GPSP'};
n = 1000;        % Signal dimension
m = 500;         % Number of measurements
rho=0.5;
rflip=0.1;
Sigma=0.1:0.2:0.5;
K=5;
% n = 5000;        % Signal dimension
% m = 1000;    
% rflip = 0.01;    % p robability of sign flips
% Sigma=0.1:0.2:1;
% rho = 0.1;
%  K  =15;
% 
len = length(Sigma);
maxnumtest = 100;
nmethod=9; 
% seednum = 0;
% rand('seed',seednum);   % fix seed 
% randn('seed', seednum); % fix seed 
total_time = zeros(maxnumtest,nmethod,len);     
total_err = zeros(maxnumtest,nmethod,len);
total_supp= zeros(maxnumtest,nmethod,len);
total_hd=zeros(maxnumtest,nmethod,len);
total_he=zeros(maxnumtest,nmethod,len);
for k = 1:len
    k
    sigma=Sigma(k)
 
    for ii = 1 : maxnumtest
        ii   
        % --------- Generate data ---------------------
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
           SIGMA(i,j) = rho^(abs(i-j));
        end
    end
    Mu = zeros(1,n);
    Psi =  mvnrnd(Mu,SIGMA,m);
    ye = Psi*xtrue;
    noise = sigma*randn(m,1);
    y = sign( ye + noise);
    nflip = floor(rflip*m);
    indxflip = randperm(m);
    indxflip = indxflip(1:nflip);
    y(indxflip) = -y(indxflip);
    %--------------------------------------------------------
 
  
   
   ff='BIHT';
   if ismember(ff,FIELDS) 
    fprintf('------BIHT ------\n')
    tic
    x = BIHT_1(y, Psi, K);
    bihttime=toc;
    mm=1;    
    total_time(ii,mm,k) = bihttime;
    total_err(ii,mm,k) =  norm(xtrue - x);
    total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
    total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,mm,k) =  1;
        end
    end
   end
  
   ff='AOP';
   if ismember(ff,FIELDS) 
    fprintf('------AOP ------\n')
    L  = nflip;  % using the real number of wrong labels.
    alpha   = 1;
    tic
    x = BIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha);
    x = x/norm(x);
    aoptime=toc;
    mm=2;
    total_time(ii,mm,k) = aoptime;
    total_err(ii,mm,k) =  norm(xtrue - x);
    total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
    total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,mm,k) = 1;
        end
    end
   end  
    %-------------------------------------------------------
    ff='LinProj';
    if ismember(ff,FIELDS) 
    fprintf('------ LinProj ------\n')
       tic,
     tx = Psi'*y/m;
     x = LinProj(tx,sqrt(K),1,.5,100,50,1e-6);
     lptime=toc;
     mm=3;
     total_time(ii,mm,k) =  lptime;
     total_err(ii,mm,k) =  norm(xtrue - x); 
     total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
    total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
     esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,mm,k) = 1;
        end
    end
    end
    %---------------------------------------------------
  ff='PBAOP'; 
   if ismember(ff,FIELDS) 
    fprintf('------ PBAOP ------\n')
     L = nflip;      % using  number of sign filips.
    alpha   = 1;
    tau = 0.05;
    tic
    x = PIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha, tau);
    x = x/norm(x);
    pbtime=toc;
    mm=4;
   total_time(ii,mm,k) = pbtime;
    total_err(ii,mm,k) =  norm(xtrue - x);
    total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
    total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,mm,k) = 1;
        end
    end
   end
   
    ff='PDASC';
   if ismember(ff,FIELDS) 
   fprintf('------L1-LS(PDASC) ------\n')
        tic
           [x,lam,ithist] = pdasc(Psi,Psi',y);
           esupp = find(x);
           Psiesupp = Psi(:,esupp); 
           % debias
           x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
           pdastime=toc;
           mm=5;
           total_time(ii,mm,k) =  pdastime;          
           total_err(ii,mm,k) =  norm(xtrue - x);
           total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
           total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
           rs = setdiff(esupp,supp);
           rse = setdiff(supp,esupp);
           if isempty(rs) && isempty(rse)
                total_supp(ii,mm,k) = 1;
           end
   end
%--------------------------------------------------------
  ff='WPDASC';
   if ismember(ff,FIELDS) 
   fprintf('------ WPDASC------\n')
         tic
        [x,lam,ithist] = wpdasc(Psi,Psi',y);
        esupp = find(x);
         Psiesupp = Psi(:,esupp); 
         % debias
         x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
         x = x/norm(x);
         pdasbtime=toc;
         mm=6;
         total_time(ii,mm,k) =  pdasbtime;          
         total_err(ii,mm,k) =  norm(xtrue - x);
         total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
         total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
         esupp = find(x);
       if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
          if isempty(rs)
             total_supp(ii,mm,k) = 1;
          end
        end
   end
 
   
   % ------------------------------------------------
    ff='MCP';
      if ismember(ff,FIELDS) 
       fprintf('------ MCP------\n')
      
       tic
        x=mcp_1bit(y,Psi,xtrue,K); 
      mcptime=toc;
        mm=7;
         total_time(ii,mm,k) =  mcptime;          
         total_err(ii,mm,k) =  norm(xtrue - x);
         total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
         total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
         esupp = find(x);
         if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
             if isempty(rs)
             total_supp(ii,mm,k) = 1;
             end
         end
      end
      
  
      ff='GPSP';
      if ismember(ff,FIELDS) 
       fprintf('------GPSP------\n')
       c      = ceil(0.01*m); 
        tic
        out      = GPSP(Psi,y,K,c); 
        x=out.x;
             %   x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
       gpsptime=toc;
        mm=8;
         total_time(ii,mm,k) =  gpsptime;          
         total_err(ii,mm,k)=  norm(xtrue - x);
         total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
         total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
         esupp = find(x);
         if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
             if isempty(rs)
             total_supp(ii,mm,k) = 1;
             end
         end
      end    
         
        % ------------------------------------------------
    ff='MCPWP';
      if ismember(ff,FIELDS) 
       fprintf('------ MCPWP------\n')
             tic
       x=mcp_pdasc(Psi,Psi',y); 
        esupp = find(x);
           Psiesupp = Psi(:,esupp); 
           % debias
           x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
       mcpwptime=toc;
        mm=9;
         total_time(ii,mm,k) =  mcpwptime;          
         total_err(ii,mm,k)=  norm(xtrue - x);
         total_hd(ii,mm,k)=nnz(sign(Psi*x-y))/m;
         total_he(ii,mm,k)=nnz(sign(Psi*x-ye))/m;
         esupp = find(x);
         if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
             if isempty(rs)
             total_supp(ii,mm,k) = 1;
             end
         end
         
      end
    end
end
averagetime = mean(total_time);
averageerror = mean(total_err);
proboracle = mean(total_supp);



 err1 = total_err(:,:,1);
 err2 = total_err(:,:,2);
 err3 = total_err(:,:,3);
 snr1=-20*log10(err1);
 snr2=-20*log10(err2);
 snr3=-20*log10(err3);
% algs = {'A1','A2','A3','A4','A5','A6','A7','A8','A9'};
 algs = {'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','GPSP','MCPWP'};
figure(1)
boxplot(snr1,'Labels',algs)
h=xlabel('method');
set(h,'Interpreter','latex','fontsize',7)
h=ylabel('SNR');
figure(2)
boxplot(snr2,'Labels',algs)
h=xlabel('method');
set(h,'Interpreter','latex','fontsize',7)
h=ylabel('SNR');
figure(3)
boxplot(snr3,'Labels',algs)
set(h,'Interpreter','latex','fontsize',7)
h=xlabel('method');
h=ylabel('SNR');


% hd1 = total_hd(:,:,1);
% hd2 = total_hd(:,:,2);
%  hd3 = total_hd(:,:,3);
%  algs = {'A1','A2','A3','A4','A5','A6','A7','A8','A9'};
% % algs = {'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP','GPSP'};
% figure(1)
% boxplot(hd1,'Labels',algs)
% h=xlabel('method');
% h=ylabel('HD');
% figure(2)
% boxplot(hd2,'Labels',algs)
% h=xlabel('method');
% h=ylabel('HD');
% figure(3)
% boxplot(hd3,'Labels',algs)
% h=xlabel('method');
% h=ylabel('HD');


 %proboracle(:,:,9)
%  time=[averagetime(:,:,1);averagetime(:,:,2);averagetime(:,:,3)];
% xlswrite('C:\installSoftware\matlab_code\sigmaCom1.xls',[prob,err,time], 'sheet1');%s     
% 
% plot(Sigma,proboracle(:,1,:);proboracle(:,2,:);proboracle(:,3,:);proboracle(:,4,:);proboracle(:,5,:);proboracle(:,6,:);proboracle(:,7,:);proboracle(:,8,:)]';%proboracle(:,:,9)
% figure(1)
% plot(Sigma,proboracle(:,1,:),'kd-',Sigma,proboracle(:,2,:),'b*:',Sigma,proboracle(:,3,:),'ro--',Sigma,proboracle(:,4,:),'x-.',proboracle(:,5,:),'gs:',Sigma,proboracle(:,6,:),'m+-',Sigma,proboracle(:,7,:),'bp-',Sigma,proboracle(:,8,:),'ch-.',Sigma,proboracle(:,9,:),'rd-','LineWidth',2)%,prob(:,9),
%  h = xlabel('$\sigma$');
% set(h,'Interpreter','latex','fontsize',12)
% set(gca,'xtick',0:0.2:1)
% h = ylabel('Probability');
% set(gca,'ytick',0:0.2:1.4)
% set(h,'Interpreter','latex','fontsize',12)
% h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP'},'FontSize',8);%,'$Homotopy$'
% set(h,'Interpreter','latex','fontsize',8)
% axis([0 1 0 1.6])
% 
% err = [averageerror(:,1,:);averageerror(:,2,:);averageerror(:,3,:);averageerror(:,4,:);averageerror(:,5,:);averageerror(:,6,:);averageerror(:,7,:);averageerror(:,8,:)]';%;averageerror(:,:,9)
% figure(2)
% plot(Sigma,err(:,1),'kd-',Sigma,err(:,2),'b*:',Sigma,err(:,3),'ro--',Sigma,err(:,4),'x-.',Sigma,err(:,5),'gs:',Sigma,err(:,6),'m+-',Sigma,err(:,7),'bp-',Sigma,err(:,8),'ch-.','LineWidth',2)%'rd-',Sigma,err(:,9),
% h = xlabel('$\sigma$');
% set(h,'Interpreter','latex','fontsize',12)
% set(gca,'xtick',0:0.2:1)
% h = ylabel('$\ell_2$-error');
% set(gca,'ytick',0:0.3:1.8)
% set(h,'Interpreter','latex','fontsize',12)
% h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP'},'FontSize',8);%,'$Homotopy$'
% set(h,'Interpreter','latex','fontsize',8)
% axis([0 1 0 1.8])

% h = xlabel('$\sigma$');
% set(h,'Interpreter','latex','fontsize',12)
% set(gca,'xtick',0:0.2:1)
% h = ylabel('$\ell_2$-error');
% set(gca,'ytick',0:0.2:1.4)
% set(h,'Interpreter','latex','fontsize',12)
% h = legend({'PDASC','BIHT','AOP','LinProj','PBAOP','WPDASC','Passive','MCP'},'FontSize',8);%,'$Homotopy$'
% set(h,'Interpreter','latex','fontsize',8)
% axis([0 1 0 1.5])

% a=[0  0.9900    0.2900    0.9800    0.7800    0.9700    0.5800    0.9600 0.0300;
%    0    0.8900    0.3200    0.8800    0.7800    0.9400    0.6700    0.8900  0.0400;
%    0    0.4000    0.3100    0.3500    0.5900    0.8700    0.5800    0.8200  0.05];
%  plot(Sigma,a(:,1),'kd-',Sigma,a(:,2),'b*:',Sigma,a(:,3),'ro--',Sigma,a(:,4),'x-.',Sigma,a(:,5),'gs:',Sigma,a(:,6),'m+-',Sigma,a(:,7),'bp-',Sigma,a(:,8),'ch-.',Sigma,a(:,9),'rd-','LineWidth',2)
% h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP','GPSP'},'FontSize',5)
%   axis([0 0.6 0 1.2])
% 
%   b=[0.0844    0.2746    0.0092    0.0960    0.0277    0.0870    0.0710    0.0302   1.3348;
%     0.0957    0.3049    0.0095    0.1056    0.0311    0.0890    0.0752     0.0263   1.3995; 
%     0.1013    0.2756    0.0085    0.0926    0.0234    0.0643    0.0662    0.0219    1.3099];
%  plot(Sigma,b(:,1),'kd-',Sigma,b(:,2),'b*:',Sigma,b(:,3),'ro--',Sigma,b(:,4),'x-.',Sigma,b(:,5),'gs:',Sigma,b(:,6),'m+-',Sigma,b(:,7),'bp-',Sigma,b(:,8),'ch-.',Sigma,b(:,9),'rd-','LineWidth',2)
% h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','MCPWP','GPSP'},'FontSize',3)
%   axis([0 0.6 0 1.6])