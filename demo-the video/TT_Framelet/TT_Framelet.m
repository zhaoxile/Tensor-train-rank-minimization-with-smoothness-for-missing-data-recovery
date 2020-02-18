function [M9, Out_TT_Framelet] = TT_Framelet( data, known, Nway, opts )
    N=numel(Nway);
    Orig = opts.X0; lambda1 = opts.lambda1; lambda2 = opts.lambda2;
    alpha = opts.alpha; th = opts.th; maxit_out = opts.maxit_out; tol = opts.tol; rho = opts.rho;
    %% Initialization
    R = 256; C = 256; I1 = 4; J1 = 4;
    M9 = initialization_M(Nway,known,data);  
    
%     [~,ranktube] = SVD_MPS_Rank_Estimation(Orig,th);   % Initial TT ranks
%     [X,Y] = initialMatrix(N,Nway,ranktube);
    
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
  
    dimL = zeros(1,N-1);
    dimR = zeros(1,N-1);
    IL = 1;
    for k = 1:N-1
        dimL(k) = IL*Nway(k);
        dimR(k) = prod(Nway)/dimL(k);
        IL = dimL(k);
    end 
    
    X0 = X; Y0 = Y;  M0 = M9;  
    X = cell(1,N-1); Y = cell(1,N-1);
    
    N = length(Nway);
    relerr = [];                       
    close all;

    %% Start Time measure
    t0 = tic;
    for k = 1 : maxit_out
        Mlast9 = M9;
        %% update (X,Y)
        for n = 1:N-1
            M_Temp = reshape(M9,[dimL(n) dimR(n)]);
            X{n}   = (alpha(n)*M_Temp*Y0{n}' + rho*X0{n})*pinv( alpha(n)*Y0{n}*Y0{n}' + rho*eye(size(Y0{n}*Y0{n}')));
            Y{n}   = pinv(alpha(n)*X{n}'*X{n} + rho*eye(size(X{n}'*X{n})))*(alpha(n)*X{n}'*M_Temp +rho*Y0{n});
        end
        
        %% update M by ADMM
        [M9,CZ_l1] = update_Framelet_M(data,known,Nway,X,Y,M0,opts);
        
        %% Calculate relative error
        relerr(k) = abs(norm(M9(:)-Mlast9(:)) / norm(Mlast9(:)));
    
        X0=X; Y0=Y; M0=M9;
        %% check stopping criterion
        if relerr(k)<=tol
            break
        end
    end

    %% Stop Time measure
    time = toc(t0);
    Out_TT_Framelet.time = time;
    Out_TT_Framelet.relerr = relerr;
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
        %uniform distribution on the unit sphere
        X0{k} = bsxfun(@rdivide,X0{k},sqrt(sum(X0{k}.^2,1)));
        Y0{k} = bsxfun(@rdivide,Y0{k},sqrt(sum(Y0{k}.^2,2)));
        IL = dimL(k);
    end   
end