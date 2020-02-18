function [M9,CZ_l1] = update_Framelet_M(data,known,Nway,X,Y, M0,opts)
    alpha = opts.alpha; beta0 = opts.beta0; lambda = opts.lambda; maxit_in = opts.maxit_in; sigma = opts.sigma; rho = opts.rho;
    frame=opts.frame; Level =opts.Level; wLevel=opts.wLevel;max_sigma=opts.max_sigma;max_beta=opts.max_beta;
    
    n1=256; n2 = 256; I1 = 2; J1 = 2;
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
    X_transition=Unfold(CastKet2Image22(M0,n1,n2,I1,J1),X_way,1);

    Z=FraDecMultiLevel(X_transition,D,Level);
    E=FraDecMultiLevel(zeros(size(X_transition)),D,Level);%(E/sigma)  
    
    for r = 1:maxit_in
          %% update M
            C1 = CoeffOper('*',E,1/sigma);
            C = CoeffOper('+',Z,C1);
            tempC=Unfold(Fold(FraRecMultiLevel(C,R,Level),X_way,1),X_way,3);
            M2 = sigma*tempC+rho*Unfold(CastKet2Image22(M0,n1,n2,I1,J1),X_way,3);
            
            for n = 1:N-1
                M2 = M2+Unfold(CastKet2Image22(beta(n)*A{n}+CC{n},n1,n2,I1,J1),X_way,3);
            end
            M3 = Fold(M2,X_way,3);
            
            for i = 1:X_way(end)
                M3(:,:,i) = M3(:,:,i)*pinv( ( sum(beta)+sigma+rho )*eye(n1,n2) );
            end                   
            M3 = min( ( max( M3, 0 ) ), 255 );
            M9 = CastImageAsKet22( M3, Nway, I1, J1);
            M9(known) = data;
           %% update A
            for n = 1:N-1
                A{n} = alpha(n)*reshape(X{n}*Y{n},Nway)+beta(n)*M9-CC{n};
                A{n} = A{n}/(alpha(n)+beta(n));
            end
           %% update Z
             M3 = CastKet2Image22(M9,n1,n2,I1,J1);
             X_transition=Unfold(M3,X_way,1);
             WMT=FraDecMultiLevel(X_transition,D,Level);
             Thresh=lambda/sigma;
             CZ_l1 = 0;
           
             for ki=1:Level
                 for ji=1:nD-1
                     for jj=1:nD-1
                         Z{ki}{ji,jj}=wthresh(WMT{ki}{ji,jj}-E{ki}{ji,jj}/sigma,'s',Thresh);
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
         %% update CC,E
          for n =1:N-1
               CC{n} = CC{n}+beta(n)*(A{n}-M9 );
          end
          for ki=1:Level
              for ji=1:nD-1
                   for jj=1:nD-1
                      if ((ji~=1)||(jj~=1))||(ki==Level)
                         deltab=Z{ki}{ji,jj}-WMT{ki}{ji,jj};
                         E{ki}{ji,jj}=E{ki}{ji,jj}+sigma*deltab;
                      end
                   end
              end
          end 
    end
end