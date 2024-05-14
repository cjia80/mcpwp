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

Sigma=0.3;
K=5;
% n = 5000;        % Signal dimension
% m = 1000;    
% rflip = 0.01;    % p robability of sign flips
% Sigma=0.1:0.2:1;
% rho = 0.1;
%  K  =15;
flip=0.05:0.05:0.15;
% 
len = length(flip);
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
    sigma=Sigma;    
 rflip=flip(k);
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
%   err4 = total_err(:,:,4);
%  err5 = total_err(:,:,5);
 snr1=-20*log10(err1);
 snr2=-20*log10(err2);
 snr3=-20*log10(err3); 
%  snr4=-20*log10(err4);
%  snr5=-20*log10(err5);
% algs = {'A1','A2','A3','A4','A5','A6','A7','A8','A9'};
 algs = {'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','GPSP','MCPWP'};
figure(1)
boxplot(snr1,'Labels',algs)
h=xlabel('method');
set(h,'Interpreter','latex','fontsize',9)
h=ylabel('SNR');

figure(2)
boxplot(snr2,'Labels',algs)
h=xlabel('method');
set(h,'Interpreter','latex','fontsize',9)
h=ylabel('SNR');
figure(3)
boxplot(snr3,'Labels',algs)
set(h,'Interpreter','latex','fontsize',9)
h=xlabel('method');
h=ylabel('SNR');


a=[1.4222    2.6076    0.0397    1.3378    0.4143    1.0127    0.3066    7.9704  0.3949;
   2.1235    3.4696    0.0403    1.7129    0.4601    1.0015    0.3205    9.1824    0.4084;
   0.4883    0.9870    0.0211    0.4288    0.1061    0.3611    0.1739    3.8899  0.1310];
   K=0.05:0.05:0.15;
figure(4) 
axis([0.05 0.15 0  10])
hold on 
 plot(K,a(:,1),'kd-',K,a(:,2),'b*:',K,a(:,3),'ro--',K,a(:,4),'x-.',K,a(:,5),'s-',K,a(:,6),'bp-',K,a(:,7),'ch-.',K,a(:,8),'gs:',K,a(:,9),'m+-','LineWidth',2)
 h = xlabel(' the sign flip ');
set(h,'Interpreter','latex','fontsize',8)
set(gca,'xtick',0.05:0.05:0.15)
h = ylabel('Times');
set(gca,'ytick',0:2:10)
set(h,'Interpreter','latex','fontsize',8)
h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','GPSP', 'MCPWP'},'FontSize',8);%,'$Homotopy$'
set(h,'Interpreter','latex','fontsize',8)

a=[ 0.1278    0.3727    0.0119    0.1430    0.0482    0.1355    0.0995    0.1274 0.0522;
   0.1327   0.3856    0.0117    0.1452    0.0422    0.1212    0.1011    0.1289  0.0442;
  0.1198    0.3720    0.0113    0.1428    0.0379    0.1137    0.1004    0.1313 0.0433];
   K=0.05:0.05:0.15;
figure(4)
 plot(K,a(:,1),'kd-',K,a(:,2),'b*:',K,a(:,3),'ro--',K,a(:,4),'x-.',K,a(:,5),'s-',K,a(:,6),'bp-',K,a(:,7),'ch-.',K,a(:,8),'gs:',K,a(:,9),'m+-','LineWidth',2)
 hold on 
 h = xlabel(' the sign flip ');
set(h,'Interpreter','latex','fontsize',8)
set(gca,'xtick',0.05:0.05:0.15)
h = ylabel('Times');
set(gca,'ytick',0:0.05:0.2)
set(h,'Interpreter','latex','fontsize',8)
h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','GPSP', 'MCPWP'},'FontSize',8);%,'$Homotopy$'
set(h,'Interpreter','latex','fontsize',8)
axis([0.05 0.15 0  0.2])


a=[0.0200    0.8600    0.5000    0.7800    0.9300    1.0000    0.7700    0.8900  0.99;
    0    0.8900    0.3200    0.8800    0.7800    0.9400    0.6700    0.7800  0.89;
    0    0.6800    0.1900    0.6500    0.5100    0.7300    0.5500    0.6400  0.66];
   K=0.05:0.05:0.15;
figure(5) 
hold on 
 plot(K,a(:,1),'kd-',K,a(:,2),'b*:',K,a(:,3),'ro--',K,a(:,4),'x-.',K,a(:,5),'s-',K,a(:,6),'bp-',K,a(:,7),'ch-.',K,a(:,8),'gs:',K,a(:,9),'m+-','LineWidth',2)
 h = xlabel(' the sign flip ');
set(h,'Interpreter','latex','fontsize',8)
set(gca,'xtick',0.05:0.05:0.15)
h = ylabel('PrE');
set(gca,'ytick',0:0.2:1)
set(h,'Interpreter','latex','fontsize',8)
h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','MCP','GPSP', 'MCPWP'},'FontSize',8);%,'$Homotopy$'
set(h,'Interpreter','latex','fontsize',8)
axis([0.05 0.15 0  1.4])


