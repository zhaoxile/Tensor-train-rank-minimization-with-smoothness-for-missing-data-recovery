function [M, Out_TMac_TT] = TMac_TT( data, known, Nway, opts )
    
    Orig = opts.Otrue;
    alpha = opts.alpha; th = opts.th; maxit = opts.maxit; tol = opts.tol;
    %% Initialization    
    N = length(Nway);
    M = initialization_M(Nway,known,data);  
%     save M.mat M
%     
%     [~,ranktube] = SVD_MPS_Rank_Estimation(Orig,th);   % Initial TT ranks   
%     [X,Y] = initialMatrix(N,Nway,ranktube);
% 
%     switch th
%         case 0.01
%             save X_TT_01.mat X
%             save Y_TT_01.mat Y
%         case 0.02 
%             save X_TT_02.mat X
%             save Y_TT_02.mat Y
%         case 0.03
%             save X_TT_03.mat X
%             save Y_TT_03.mat Y
%     end 

    switch th
        case 0.01
            load X_TT_01.mat
            load Y_TT_01.mat
        case 0.02
            load X_TT_02.mat
            load Y_TT_02.mat
        case 0.03
            load X_TT_03.mat
            load Y_TT_03.mat
    end
     
    R  = 256; C = 256; I1 = 4; J1 = 4;
    Out_TMac_TT = [];
    Xsq = cell(1,N-1);
    k = 1;
    relerr = [];
    relerr(1) = 1;
    
%     close all;
    %% Initialize figures
%     h1 = figure();clf;
%     subplot(1,3,1);cla;hold on;
%     subplot(1,3,2);
%     set(h1,'Position',[200 600 1500 350]);
%     Img = CastKet2Image22(M,R,C,I1,J1);
%     imagesc(uint8(Img));drawnow;
%     cla;hold on;
%     subplot(1,3,3);
%     fprintf('iter: RSE  \n');
    
    %% Start Time measure
    t0=tic;
    while relerr(k) > tol
        k = k+1;
        Mlast = M;

        %% update (X,Y)
        for n = 1:N-1
            Mn     = reshape(M,[size(X{n},1) size(Y{n},2)]);
            X{n}   = Mn*Y{n}';
            Xsq{n} = X{n}'*X{n};
            Y{n}   = pinv(Xsq{n})*X{n}'*Mn;
        end
        %% update M
        Mn = X{1}*Y{1};
        M  = alpha(1)*Mn;
        M  = reshape(M,Nway);
        for n = 2:N-1
            Mn = X{n}*Y{n};
            Mn = reshape(Mn,Nway);
            M  = M+alpha(n)*Mn;
        end
        M(known) = data;
        %% Calculate relative error
        relerr(k) = abs(norm(M(:)-Mlast(:)) / norm(Mlast(:)));
        
        %% Update figures
%         set(0,'CurrentFigure',h1);
%         subplot(1,3,1);cla;hold on;   
%         plot(relerr);
%         plot(tol*ones(1,length(relerr)));
%         grid on
%         set(gca,'YScale','log')
%         title('Relative Error')
%         ylim([(tol-5e-5) 1]);
%         subplot(1,3,2);
        Img = CastKet2Image44(M,R,C,I1,J1);
%         imagesc(uint8(Img));
%         title('TMac-TT');        
%         drawnow; 
%         % message log
        Img0=CastKet2Image44(Orig,R,C,I1,J1);
%         subplot(1,3,3);
%         imagesc(uint8(Img0));
%         title('Original image');
%         if k == 2 || mod(k-1,200) ==0 
%               fprintf('  %d:  %f \n',k-1,RSE(Img(:),Img0(:)));
%         end
        %% check stopping criterion
        if k > maxit || (relerr(k)-relerr(k-1) > 0)
            break
        end
    end
    
    %% Stop Time measure
    time = toc(t0);
    Out_TMac_TT.time = time;
    Out_TMac_TT.relerr = relerr;
end
function [X0,Y0] = initialMatrix(N,Nway,ranktube)
    X0 = cell(1,N-1);Y0 = cell(1,N-1);
    dimL = zeros(1,N-1);
    dimR = zeros(1,N-1);
    IL = 1;
    for k = 1:N-1
        dimL(k) = IL*Nway(k);
        dimR(k) = prod(Nway)/dimL(k);
        %
        X0{k} = randn(dimL(k),ranktube(k));
        Y0{k} = randn(ranktube(k),dimR(k));
        % uniform distribution on the unit sphere
        X0{k} = bsxfun(@rdivide,X0{k},sqrt(sum(X0{k}.^2,1)));
        Y0{k} = bsxfun(@rdivide,Y0{k},sqrt(sum(Y0{k}.^2,2)));
        %
        IL = dimL(k);
    end   
end