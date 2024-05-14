% The computational time and recovery accuracy and 
% support recovery probability for  BIHT, AOP, LinProj
%,PBAOP, L1-LS(PDASC) and LQ_LS(WPDASC) and  are compared. 
clc
clear all 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
n = 10000;          % Signal dimension
m = 2000;           % Number of measurements
K = 30;             % signal sparsity
% rflip = 0.03;     % probability of sign flips
% sigma = 0.2;      % noise level
 rflip = 0.05;       % probability of sign flips
sigma = 0.3;        % noise level 
% rflip = 0.1;       % probability of sign flips
% sigma = 0.5;        % noise level 
rho=0.3;
% rflip = 0.1;     % probability of sign flips
% sigma = 0.3;      % noise level 
total_err = zeros(9,1);
total_time = zeros(9,1);
total_supp = zeros(9,1);
 total_hd=zeros(9,1);
 total_he=zeros(9,1);
seednum = 0;
rand('seed',seednum);   % fix seed 
randn('seed', seednum); % fix seed 
for ii = 1 : 1
    ii   
    % --------- Generate data ---------------------
     rand('seed',seednum);   % fix seed 
        randn('seed', seednum); % fix seed 
        xtrue = zeros(n,1);
        rp = randperm(n);
        xtrue(rp(1:K)) = sign(randn(K,1));
        xtrue = xtrue/norm(xtrue);
        supp = find(xtrue);
        SIGMA = zeros(n,n);
        for kk = 1:n
            for jj = 1:n
               SIGMA(kk,jj) = rho^(abs(kk-jj));
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
    % ------------------------------------------------
    fprintf('------ BIHT  ----------\n')
    tic
    x = BIHT_1(y, Psi, K);
    total_time(1) = total_time(1) + toc;
    total_err(1) = total_err(1) + norm(xtrue - x);
     total_hd(1)=total_hd(1)+nnz(sign(Psi*x-y))/m;
       total_he(1)=total_he(1)+nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(1) = total_supp(1) + 1;
        end
    end
    fprintf('------ AOP  ------\n')
    L  = nflip;  % using the real number of wrong labels.
    alpha   = 1;
    tic
    x = BIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha);
    x = x/norm(x);
    total_time(2) = total_time(2) + toc;
    total_err(2) = total_err(2) + norm(xtrue - x);
     total_hd(2)=total_hd(2)+nnz(sign(Psi*x-y))/m;
       total_he(2)=total_he(2)+nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(2) = total_supp(2) + 1;
        end
    end
    %-------------------------------------------------------
    fprintf('------ LinProj ----------\n')
     tic,
     tx = Psi'*y/m;
     x = LinProj(tx,sqrt(K),1,.5,100,50,1e-6);
     total_time(3) = total_time(3) + toc;
     total_err(3) = total_err(3) + norm(xtrue - x);   
      total_hd(3)=total_hd(3)+nnz(sign(Psi*x-y))/m;
       total_he(3)=total_he(3)+nnz(sign(Psi*x-ye))/m;
     esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(3) = total_supp(3) + 1;
        end
    end
    %---------------------------------------------------
    fprintf('------ PBAOP  ------\n')
    L = nflip;      % using  number of sign filips.
    alpha   = 1;
    tau = 0.05;
    tic
    x = PIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha, tau);
    x = x/norm(x);
    total_time(4) = total_time(4) + toc;
    total_err(4) = total_err(4) + norm(xtrue - x);
     total_hd(4)=total_hd(4)+nnz(sign(Psi*x-y))/m;
       total_he(4)=total_he(4)+nnz(sign(Psi*x-ye))/m;
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(4) = total_supp(4) + 1;
        end
    end
%--------------------------------------------------------
    fprintf('------ L1-LS(PDASC) ------\n')
        tic
       [x,lam,ithist] = pdasc(Psi,Psi',y);
       % debias
       esupp = find(x);
       Psiesupp = Psi(:,esupp); 
       x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
       x = x/norm(x);
       total_time(5) = total_time(5) + toc;
       total_err(5) = total_err(5) + norm(xtrue - x);
        total_hd(5)=total_hd(5)+nnz(sign(Psi*x-y))/m;
       total_he(5)=total_he(5)+nnz(sign(Psi*x-ye))/m;
       esupp = find(x);
        if length(esupp) == length(supp)
            rs = setdiff(esupp,supp);
            if isempty(rs)
                 total_supp(5) = total_supp(5) + 1;
            end
        end
        
 
         fprintf('------ GPSP------\n')
    tic
       k          = ceil(0.01*m);
      out        = GPSP(Psi,y,K,k);
      x=out.x;
             % debias
       esupp = find(x);
%        Psiesupp = Psi(:,esupp); 
%        x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
       x = x/norm(x);
       total_time(6) = total_time(6) + toc;
       total_err(6) = total_err(6) + norm(xtrue - x);   
        total_hd(6)=total_hd(6)+nnz(sign(Psi*x-y))/m;
       total_he(6)=total_he(6)+nnz(sign(Psi*x-ye))/m;
        if length(esupp) == length(supp)
            rs = setdiff(esupp,supp);
            if isempty(rs)
                 total_supp(6) = total_supp(6) + 1;
            end
        end
        
          
        %-------------------------------------------------------
          fprintf('------ LQ-LS(WPDASC) ------\n')
    tic
       [x,lam,ithist] = wpdasc(Psi,Psi',y);
       % debias
       esupp = find(x);
       Psiesupp = Psi(:,esupp); 
       x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
       x = x/norm(x);
       total_time(7) = total_time(7) + toc;
       total_err(7) = total_err(7) + norm(xtrue - x);
        total_hd(7)=total_hd(7)+nnz(sign(Psi*x-y))/m;
       total_he(7)=total_he(7)+nnz(sign(Psi*x-ye))/m;
        if length(esupp) == length(supp)
            rs = setdiff(esupp,supp);
            if isempty(rs)
                 total_supp(7) = total_supp(7) + 1;
            end
        end
         fprintf('------ MCP ------\n')
    tic
       x=mcp_1bit(y,Psi,xtrue,K);
       % debias
       esupp = find(x);
%        Psiesupp = Psi(:,esupp); 
%        x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
       x = x/norm(x);
       total_time(8)= total_time(8) + toc;
       total_err(8) = total_err(8) + norm(xtrue - x);  
        total_hd(8)=total_hd(8)+nnz(sign(Psi*x-y))/m;
       total_he(8)=total_he(8)+nnz(sign(Psi*x-ye))/m;
        if length(esupp) == length(supp)
            rs = setdiff(esupp,supp);
            if isempty(rs)
                 total_supp(8) = total_supp(8) + 1;
            end
        end 
      
        %-------------------------------------------------------
          fprintf('------ MCP-LS(MCPWP) ------\n')
    tic
       [x,lam,ithist] = mcp_pdasc(Psi,Psi',y);
       % debias
       esupp = find(x);
       %Psiesupp = Psi(:,esupp); 
      % x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
       x = x/norm(x);
       total_time(9) = total_time(9) + toc;
       total_err(9) = total_err(9) + norm(xtrue - x);  
       err=norm(xtrue - x)
        total_hd(9)=total_hd(9)+nnz(sign(Psi*x-y))/m;
       total_he(9)=total_he(9)+nnz(sign(Psi*x-ye))/m;
        if length(esupp) == length(supp)
            rs = setdiff(esupp,supp);
            if isempty(rs)
                 total_supp(9) = total_supp(9) + 1;
            end
        end
end

% show results 
fprintf('------ average  error ----------\n')
err=total_err/ii
fprintf('------ average  snr ----------\n')
-20*log10(total_err/ii)
fprintf('------ average   time ----------\n')
total_time/ii
fprintf('------  support recovery probability ----------\n')
total_supp/ii
 total_hd/ii
 total_he/ii
