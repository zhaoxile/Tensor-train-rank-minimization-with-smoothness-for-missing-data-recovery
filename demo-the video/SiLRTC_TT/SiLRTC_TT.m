function [M, Out_Si_TT] = SiLRTC_TT( data, known, Nway, opts )

    Orig  = opts.Otrue;
    alpha = opts.alpha; beta = opts.beta; maxit = opts.maxit; tol = opts.tol;
   %% Initialization 
    N = length(Nway);
    M = initialization_M(Nway,known,data);

    
    Out_Si_TT = [];
   
    
    dimL = zeros(1,N-1);
    dimR = zeros(1,N-1);
    IL = 1;
    for m = 1:N-1
        dimL(m) = IL*Nway(m);
        dimR(m) = prod(Nway)/dimL(m);
        IL = dimL(m);
    end 
    
    X = cell(1,N-1);
    k = 1;
    relerr = [];
    relerr(1) = 1;
 
    close all;
    % Initialize figures
    R=256;C = 256;I1 = 4;J1 = 4;
%     h1 = figure();clf;
%     subplot(1,3,1);cla;hold on;
%     subplot(1,3,2);
%     set(h1,'Position',[200 600 1500 350]);
%     Img = CastKet2Image22(M,R,C,I1,J1);
%     Img = reshape(Img,[R C 3]);
%     imagesc(uint8(Img));drawnow;
%     cla;hold on;
%     subplot(1,3,3);
%     fprintf('iter: RSE  \n');
    % Start Time measure
    t0=tic;
    while relerr(k) > tol
        k = k+1;
        Mlast = M;

       %% update X
        for n = 1:N-1
            Xn   = reshape(M, [dimL(n) dimR(n)]);
            X{n} = SVT( Xn, alpha(n)/beta(n));
        end
        
        %% update M
        M = beta(1)*X{1};
        M = reshape(M, Nway);
        for n = 2:N-1
            Mn = reshape(X{n}, Nway);
            M  = M+beta(n)*Mn;
        end
        M = M./(sum(beta));
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
%         title('SiLRTC-TT');        
%         drawnow; 
%         % message log
        Img0=CastKet2Image44(Orig,R,C,I1,J1);
%         subplot(1,3,3);
%         imagesc(uint8(Img0));
%         title('Original image');
%         if k == 2 || mod(k-1,50) ==0 
%             fprintf('  %d:  %f \n',k-1,RSE(Img(:),Img0(:)));
%         end 
        %% check stopping criterion
        if k > maxit || (relerr(k)-relerr(k-1) > 0)
            break
        end
    end
    % Stop Time measure
    time = toc(t0);
    Out_Si_TT.time = time;
    Out_Si_TT.relerr = relerr;
end