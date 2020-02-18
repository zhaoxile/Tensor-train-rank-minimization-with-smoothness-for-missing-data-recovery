function X_out=Framelet_X3(A,Y,X_k,rho,opts2)

if isfield(opts2,'mu');         mu = opts2.mu;                   else  mu = 10;                   end
if isfield(opts2,'beta');   beta= opts2.beta;                     else  beta = 1000;               end
if isfield(opts2,'F_tol');       tol = opts2.F_tol;                  else  tol = 1e-5;                end
if isfield(opts2,'frame');      frame = opts2.frame;           else  frame = 1;                  end
if isfield(opts2,'Level');        Level = opts2.Level;               else  Level= 1;                    end
if isfield(opts2,'wLevel');      wLevel = opts2.wLevel;          else  wLevel= -1;                 end
if isfield(opts2,'F_it');         maxit = opts2.F_it;                else  maxit = 15;                  end
if isfield(opts2,'x_size');      x_size =opts2.x_size;            else  x_size = sqrt(size(Y,2));end

[D,R]=GenerateFrameletFilter(frame);
nD=length(D);
m=x_size(1);
n=x_size(2);
normg=norm(Y,'fro');
e=size(A,2);
X_way=[m,n,e];
X_transition=zeros(m,n*e);
X_out=zeros(e,m*n);

Z=FraDecMultiLevel(X_transition,D,Level);

Theta=Z;
for nstep=1:maxit
    
    %% solve X
    temp=mu*A'*Y;
    C=CoeffOper('-',Z,Theta);
    tempC=Unfold(Fold(FraRecMultiLevel(C,R,Level),X_way,1),X_way,3);
    X_out=pinv(mu*A'*A+(rho+beta)*eye(size(A'*A)))*(rho*X_k+temp+beta*tempC); %#ok<MHERM>
    
    %% solve Z
    X_transition=Unfold(Fold(X_out,X_way,3),X_way,1);
    C=FraDecMultiLevel(X_transition,D,Level);
    Thresh=1/beta;
    for ki=1:Level
        for ji=1:nD-1
            for jj=1:nD-1
                Z{ki}{ji,jj}=wthresh(C{ki}{ji,jj}+Theta{ki}{ji,jj},'s',Thresh);
            end
        end
        if wLevel<=0
            Thresh=Thresh*norm(D{1});
        else
            Thresh=Thresh*wLevel;
        end
    end
    %% update Theta
    error=0;
    for ki=1:Level
        for ji=1:nD-1
            for jj=1:nD-1
                if ((ji~=1)||(jj~=1))||(ki==Level)
                    deltab=C{ki}{ji,jj}-Z{ki}{ji,jj};
                    error=error+norm(deltab,'fro')^2;
                    Theta{ki}{ji,jj}=Theta{ki}{ji,jj}+deltab;
                end
            end
        end
    end
    error=sqrt(error)/normg;
    if error<tol
        break;
    end
end