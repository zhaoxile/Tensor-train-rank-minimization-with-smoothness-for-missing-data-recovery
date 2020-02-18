clc
clear all
close all
rand('seed',213412);
addpath(genpath(cd));

EN_HaLRTC   = 0;
EN_tSVD     = 0;
EN_SiLRTCTT = 0;
EN_TMacTT   = 0;
EN_TT_Framelet = 1;
methodname  = { 'HaLRTC', 'tSVD', 'SiLRTC-TT', 'TMac-TT', 'TT-Framelet' };

load bus_256_256_30.mat
X0   = bus(:,:,:);
[R C B] = size(X0); 
for i = 1:B
    X = X0(:,:,i);
    X0(:,:,i)=X./max(max(X(:)))*255;
end

name = {'bus'};
Nway = [16 16 16 16 B];     % 9th-order dimensions for KA 
I1 = 4; J1 = 4; 

for SR = [0.1]
    
%% Sampling   
sizeData = size(X0);
Y = zeros(sizeData);
Mask  = zeros(sizeData);
Index = find(rand(prod(sizeData),1)<SR);
Y(Index) = X0(Index);
Mask(Index) = 1;

Otrue  = CastImageAsKet44( X0, Nway, I1 ,J1 );
Oknown = CastImageAsKet44( Mask, Nway, I1, J1 );
Oknown = find( Oknown==1 );
Okn    = Otrue( Oknown );

Omiss = zeros( Nway );
Omiss( Oknown ) = Otrue( Oknown );
Omiss = CastKet2Image44( Omiss, 256, 256, I1, J1 );

imname=[num2str(name{1}),'_miss_SR_',num2str(SR),'.mat'];
save(imname,'X0','Omiss', 'Oknown','Okn');
%% HaLRTC
j = 1;
if EN_HaLRTC
    %%%%
    fprintf('\n');
    disp(['performing ',methodname{j}, ' ... ']);
  
    opts = [];
    alpha      = ones(ndims(X0),1);
    opts.alpha = alpha/sum(alpha); opts.tol = 1e-4; opts.maxit = 500; opts.rho = 1.1; opts.max_beta = 1e10;

    for beta = [5*1e-4]
        opts.beta = beta;
        tic;
        X = LRTC_HaLRTC( Omiss, Index, opts );
        time=toc;
        PSNRvector=zeros(1,B);
        for i=1:1:B
            PSNRvector(i)=psnr3(X0(:,:,i)/255,X(:,:,i)/255);
        end
        psnr = mean(PSNRvector);
               
        SSIMvector=zeros(1,B);
        for i=1:1:B
            SSIMvector(i)=ssim3(X0(:,:,i),X(:,:,i));
        end
         ssim  = mean(SSIMvector);
 
        display(sprintf('psnr=%.2f,ssim=%.4f,beta=%.5f', psnr, ssim, opts.beta))
        display(sprintf('=================================='))
        
        imname=[num2str(name{1}),'_SR_',num2str(SR),'_result_HaLRTC_psnr_',num2str(psnr,'%.2f'),'_ssim_',num2str(ssim,'%.4f'),'_beta_',num2str(opts.beta,'%.5f'),'.mat'];
        save(imname,'X','PSNRvector','SSIMvector','time');
    end
end
%% tSVD
j = j+1;
if EN_tSVD
    %%%%
    fprintf('\n');
    disp(['performing ',methodname{j}, ' ... ']);
    
    alpha = 1; maxItr  = 1000; myNorm = 'tSVD_1'; 
    A     = diag(sparse(double(Mask(:)))); 
    b     = A * X0(:);
    beta1 = [5*1e-6];
    
    for k = 1:length(beta1)
        beta = beta1(k);
        tic;
        X =  tensor_cpl_admm( A, b, beta, alpha, sizeData, maxItr, myNorm, 0 );
        X = reshape(X,sizeData);
        time=toc;
        PSNRvector=zeros(1,B);
        for i=1:1:B
            PSNRvector(i)=psnr3(X0(:,:,i)/255,X(:,:,i)/255);
        end
        psnr = mean(PSNRvector);   
        
        SSIMvector=zeros(1,B);
        for i=1:1:B
            SSIMvector(i)=ssim3(X0(:,:,i),X(:,:,i));
        end
         ssim  = mean(SSIMvector);
        
        display(sprintf('psnr=%.2f,ssim=%.4f,beta=%.6f', psnr, ssim, beta))
        display(sprintf('=================================='))
        
        imname=[num2str(name{1}),'_SR_',num2str(SR),'_result_tSVD_psnr_',num2str(psnr,'%.2f'),'_ssim_',num2str(ssim,'%.4f'),'_beta_',num2str(beta,'%.6f'),'.mat'];
        save(imname,'X','PSNRvector','SSIMvector','time');
    end
end
%% SiLRTC-TT
j = j + 1;
if EN_SiLRTCTT
    %%%%
    fprintf('\n');
    disp(['performing ',methodname{j}, ' ... ']);

    opts=[]; 
    opts.alpha = weightTC(Nway); opts.tol = 1e-4; opts.maxit = 1000; opts.Otrue = Otrue;
    %%%%
    r = [0.01];
    for k = 1:length(r)
        opts.beta = r(k)*opts.alpha;
        tic;
        [X, Out] = SiLRTC_TT( Okn, Oknown, Nway, opts );
         time=toc;
        X = CastKet2Image44(X,256,256,I1,J1);
        
        PSNRvector=zeros(1,B);
        for i=1:1:B
            PSNRvector(i)=psnr3(X0(:,:,i)/255,X(:,:,i)/255);
        end
        psnr = mean(PSNRvector);
        
        SSIMvector=zeros(1,B);
        for i=1:1:B
            SSIMvector(i)=ssim3(X0(:,:,i),X(:,:,i));
        end
         ssim  = mean(SSIMvector);
  
        display(sprintf('psnr=%.2f,ssim=%.4f,r=%.2f',psnr, ssim, r(k)))
        display(sprintf('=================================='))
        
        imname=[num2str(name{1}),'_SR_',num2str(SR),'_result_Si_TT_psnr_',num2str(psnr,'%.2f'),'_ssim_',num2str(ssim,'%.4f'),'_r_',num2str(r(k)),'.mat'];
        save(imname,'X','PSNRvector','SSIMvector','time');
    end
end
%% TMac-TT
j = j + 1;
if EN_TMacTT
    %%%%
    fprintf('\n');
    disp(['performing ',methodname{j}, ' ... ']);
    
    opts=[];
    opts.alpha = weightTC(Nway); opts.tol   = 1e-4; opts.maxit = 1000; opts.Otrue = Otrue;
    ps = [];
    XX = cell(1,3);
    th = [0.03];
%     opts.th = th;
    for k = 1:length(th)
        opts.th = th(k);
        tic;
        [X, Out] = TMac_TT( Okn, Oknown, Nway, opts );
        X = CastKet2Image44(X,256,256,I1,J1);
        time=toc;
        PSNRvector=zeros(1,B);
        for i=1:1:B
            PSNRvector(i)=psnr3(X0(:,:,i)/255,X(:,:,i)/255);
        end
        psnr = mean(PSNRvector);
                
        SSIMvector=zeros(1,B);
        for i=1:1:B
            SSIMvector(i)=ssim3(X0(:,:,i),X(:,:,i));
        end
         ssim  = mean(SSIMvector);
        
        ps(k) = psnr;
        XX{k} = X;

        display(sprintf('psnr=%.2f,ssim=%.4f,th=%.2f',psnr, ssim, opts.th))
        display(sprintf('=================================='))
        
        imname=[num2str(name{1}),'_SR_',num2str(SR),'_result_TMac_TT_psnr_',num2str(psnr,'%.2f'),'_ssim_',num2str(ssim,'%.4f'),'_th_',num2str(opts.th),'.mat'];
        save(imname,'X','PSNRvector','SSIMvector','time');  
    end
end

%% TT-Framelet
j = j+1;
if EN_TT_Framelet
    opts=[];
    opts.alpha  = weightTC(Nway); 
    opts.X0 = X0;
    opts.tol    = 1e-5;
    opts.maxit_out  = 100;
    opts.maxit_in  = 15;   
    opts.max_sigma1 = 10; 
    opts.max_sigma2 = 10; 
    opts.max_beta = 10; 
    opts.rho    = 10^(-3);
    opts.th     = 0.03;
    opts.lambda1 = 0.5; 
    opts.lambda2 = 0.5; 
    opts.beta0  = 0.01; 
    opts.sigma1 = 0.01;      
    opts.sigma2 = 0.001;   
    opts.frame = 3;     
    opts.Level = 1;
    opts.wLevel= 0.5;
    %%%%%
    fprintf('\n');
    disp(['performing ',methodname{j}, ' ... ']);

        tic;
        [X_TT_Framelet, Out_TT_Framelet] = TT_Framelet( Okn, Oknown, Nway, opts );
        time=toc;
        X_TT_Framelet    = CastKet2Image44(X_TT_Framelet,256,256,I1,J1);
        
        PSNRvector=zeros(1,B);
        for i=1:1:B
            PSNRvector(i)=psnr3(X0(:,:,i)/255,X_TT_Framelet(:,:,i)/255);
        end
        PSNR_TT_Framelet = mean(PSNRvector);
        
        RSEvector=zeros(1,B);
        for i=1:1:B
            RSEvector(i)=RSE(X0(:,:,i),X_TT_Framelet(:,:,i));
        end
        RSE_TT_Framelet = mean(RSEvector);
        
        SSIMvector=zeros(1,B);
        for i=1:1:B
            SSIMvector(i)=ssim3(X0(:,:,i),X_TT_Framelet(:,:,i));
        end
        SSIM_TT_Framelet = mean(SSIMvector);

        display(sprintf('psnr=%.2f,ssim=%.4f,rse=%.3f,th=%.2f,lambda1=%.3f,lambda2=%.3f,beta0=%.3f,sigma1=%.3f,sigma2=%.3f',PSNR_TT_Framelet, SSIM_TT_Framelet, RSE_TT_Framelet, opts.th, opts.lambda1, opts.lambda2, opts.beta0, opts.sigma1,opts.sigma2))
        display(sprintf('=================================='))
        
        imname=[num2str(name{1}),'_SR_',num2str(SR),'_result_TT_Framelet_psnr_',num2str(PSNR_TT_Framelet,'%.2f'),'_ssim_',num2str(SSIM_TT_Framelet,'%.4f'),'_Rse_',num2str(RSE_TT_Framelet,'%.3f'),'_th_',num2str(opts.th),'_lambda1_',num2str(opts.lambda1),'_lambda2_',num2str(opts.lambda2),'_beta0_',num2str(opts.beta0),'_sigma1_',num2str(opts.sigma1),'_sigma2_',num2str(opts.sigma2),'.mat'];
        save(imname,'X_TT_Framelet','PSNRvector','RSEvector','SSIMvector','time');       
    end
    end
