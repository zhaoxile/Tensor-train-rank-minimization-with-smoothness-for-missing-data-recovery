function [M9,CZ_l1] = update_Framelet_M(data,known,Nway,X,Y, M0,opts)
    alpha = opts.alpha; beta0 = opts.beta0; lambda1 = opts.lambda1; lambda2 = opts.lambda2; sigma1 = opts.sigma1; sigma2 = opts.sigma2;
    frame=opts.frame; Level =opts.Level; wLevel=opts.wLevel; maxit_in = opts.maxit_in; rho = opts.rho;
    max_sigma1=opts.max_sigma1;max_beta=opts.max_beta;max_sigma2=opts.max_sigma2;
    
    n1=256; n2 = 256; I1 = 4; J1 = 4;
    N=numel(Nway);
    A = cell(1,N-1); CC = cell(1,N-1);
    for i = 1:N-1
        A{i} = M0;
        CC{i} = zeros(size(A{i}));
    end
    beta=beta0*alpha;
    M9 = M0;

    [D,R]=GenerateFrameletFilter(frame);
    nD=length(D);
    X_way=[n1,n2,Nway(end)];
    X_transition=Unfold(CastKet2Image44(M0,n1,n2,I1,J1),X_way,1);
    Z=FraDecMultiLevel(X_transition,D,Level);
    E=FraDecMultiLevel(zeros(size(X_transition)),D,Level);%(E/sigma)  
    P=Unfold(CastKet2Image44(M0,n1,n2,I1,J1),X_way,3);
    F=zeros(size(P));
    [conjoDx,conjoDy,Denom1,Denom2] = getC(P);
    
    for r = 1:maxit_in
           %% update M
            C1 = CoeffOper('*',E,1/sigma1);
            C=CoeffOper('+',Z,C1);
            tempC=Unfold(Fold(FraRecMultiLevel(C,R,Level),X_way,1),X_way,3);
            b1 = sigma2*conjoDy.*fft2(P) + conjoDy.*fft2(F);
            b2 = rho*fft2(Unfold(CastKet2Image44(M0,n1,n2,I1,J1),X_way,3));
            M2 = sigma1*fft2(tempC)+b2+b1;
              
            for n = 1:N-1
                M2 = M2+fft2(Unfold(CastKet2Image44(beta(n)*A{n}+CC{n},n1,n2,I1,J1),X_way,3));
            end           
            a  = sum(beta)+sigma1+sigma2*Denom2+rho;
            M2 = M2./a;
            M2 = real(ifft2(M2));
     
            M3 = Fold(M2,X_way,3);
            M3 = min( ( max( M3, 0 ) ), 255 );
            M9 = CastImageAsKet44( M3, Nway, I1, J1);
            M9(known) = data;
             
           %% update A
            for n = 1:N-1
                A{n} = alpha(n)*reshape(X{n}*Y{n},Nway)+beta(n)*M9-CC{n};
                A{n} = A{n}/(alpha(n)+beta(n));
            end
              
           %% update Z
             M3 = CastKet2Image44(M9,n1,n2,I1,J1);
             X_transition=Unfold(M3,X_way,1);
             WMT=FraDecMultiLevel(X_transition,D,Level);
             Thresh=lambda1/sigma1;
             CZ_l1 = 0;
           
             for ki=1:Level
                 for ji=1:nD-1
                     for jj=1:nD-1
                         Z{ki}{ji,jj}=wthresh(WMT{ki}{ji,jj}-E{ki}{ji,jj}./sigma1,'s',Thresh);
                         tempCZ=WMT{ki}{ji,jj};
                         CZ_l1 = CZ_l1+sum(abs(tempCZ(:)));
                     end
                 end
                 if wLevel<=0
                    Thresh=Thresh*norm(D{1});
                 else
                 Thresh=Thresh*wLevel;
                 end
             end
           %% update P
              temp2 = diff(M2,1,1);  % row diff
              [m,n] = size(M2);
              dy = [temp2; temp2(m-1,:)];   
              P = soft_shrink(dy - F/sigma2,lambda2/sigma2);      
        
          %% update CC,E,F,G
           for n =1:N-1
               CC{n} = CC{n}+beta(n)*(A{n}-M9 );
           end
           for ki=1:Level
              for ji=1:nD-1
                   for jj=1:nD-1
                      if ((ji~=1)||(jj~=1))||(ki==Level)
                         deltab=Z{ki}{ji,jj}-WMT{ki}{ji,jj};
                         E{ki}{ji,jj}=E{ki}{ji,jj}+sigma1*deltab;
                      end
                   end
              end
           end
           F = F + sigma2 * (P - dy);
           beta = min(1.1*beta, max_beta); 
           sigma1 = min(1.1*sigma1, max_sigma1);
           sigma2 = min(1.1*sigma2, max_sigma2);   
    end
end

function [conjoDx,conjoDy,Denom1,Denom2] = getC(I)
% compute fixed quantities
sizeI = size(I);
otfDx = psf2otf([1,-1],sizeI);
otfDy = psf2otf([1;-1],sizeI);
conjoDx = conj(otfDx);
conjoDy = conj(otfDy);
Denom1 = abs(otfDx).^2;
Denom2 = abs(otfDy).^2;
end
