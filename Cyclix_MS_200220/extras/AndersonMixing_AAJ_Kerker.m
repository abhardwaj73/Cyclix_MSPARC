function [x_new] = AndersonMixing_AAJ_Kerker(x1, f1, m, iter, mix, p, B, D, spin,S)
% f = residual, m = history, x=iterate, f(x)=g(x)-x, p=AJ step
% B is the matrix to-be-inverted part of preconditioner, D is the Laplacian
% Reference:
%         Pratapa, Phanisri P., Phanish Suryanarayana, and John E. Pask. 
%         "Anderson acceleration of the Jacobi iterative method: 
%         An efficient alternative to Krylov methods for large, sparse linear systems." 
%         Journal of Computational Physics 306 (2016): 43-54.

persistent DX DX1 DX2 DF DF1 DF2 Xold Fold PF0 PF0_1 PF0_2 PF20 PF20_1 PF20_2
dp = 0; tol=1e-10; maxit=100;

if spin == 1
    if iter == 1
        PF1=AAR(B,D*f1,[],tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
        PF0=PF1;
    else
        PF1=AAR(B,D*f1,PF0,tol,maxit,[],[],[],[],[]); 
        PF0=PF1;
    end
    
    if iter > 1
        k = mod(iter-2,m)+1;

        DX(:,k) = x1 - Xold;
        DF(:,k) = f1 - Fold;
        Yk = pinv(DF'*DF)*(DF'*f1);        
        f_bar = DF*Yk;
        if iter == 2
            PF2(:,k)=AAR(B,D*f_bar,[],tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
        else
            PF2(:,k)=AAR(B,D*f_bar,PF20,tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
        end
        PF20 = PF2(:,k);
        
        if rem(iter,p)==0                          % apply every p iters 
            dp = DX*Yk + mix*PF20;
        else
            dp=0;
        end    
    end
    x_new = x1 + mix*PF0 - dp;
    Xold = x1;
    Fold = f1;
elseif spin == 2
   if iter == 1
        PF1_1=AAR(B,D*f1(1:end/2),[],tol,maxit,[],[],[],[],[]);             % The 1st preconditioned term
        PF1_1 = PF1_1 + dot(S.W,f1(1:end/2))/sum(S.W);
        PF1_2=AAR(B,D*f1(end/2+1:end),[],tol,maxit,[],[],[],[],[]);
        PF1_2 = PF1_2 + dot(S.W,f1(end/2+1:end))/sum(S.W);
        PF1 = vertcat(PF1_1,PF1_2);
        PF0=PF1;
    else
        PF1_1=AAR(B,D*f1(1:end/2),PF0(1:end/2),tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
        PF1_1 = PF1_1 + dot(S.W,f1(1:end/2))/sum(S.W);
        PF1_2=AAR(B,D*f1(end/2+1:end),PF0(end/2+1:end),tol,maxit,[],[],[],[],[]);
        PF1_2 = PF1_2 + dot(S.W,f1(end/2+1:end))/sum(S.W);
        PF1 = vertcat(PF1_1,PF1_2);
        PF0=PF1;
   end
%    sum(x1(1:end/2))
%    sum(x1(end/2+1:end))
   %min(f1(1:end/2))
   %min(f1(end/2+1:end))
   %min(PF1_1)
   %min(PF1_2)
   %dot(S.W,f1(1:end/2))
   %dot(S.W,f1(end/2+1:end))
   %dot(S.W,PF1_1)
   %dot(S.W,PF1_2)
    
    if iter > 1
        k = mod(iter-2,m)+1;
        DX(:,k) = x1 - Xold;
        DF(:,k) = f1 - Fold;
        Yk = pinv(DF'*DF)*(DF'*f1);
        f_bar = DF*Yk;
        if iter == 2
            PF2_1=AAR(B,D*f_bar(1:end/2),[],tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
            PF2_1 = PF2_1 + dot(S.W,f_bar(1:end/2))/sum(S.W);
            PF2_2=AAR(B,D*f_bar(end/2+1:end),[],tol,maxit,[],[],[],[],[]);
            PF2_2 = PF2_2 + dot(S.W,f_bar(end/2+1:end))/sum(S.W);
            PF2(:,k) = vertcat(PF2_1,PF2_2);            
        else
            PF2_1=AAR(B,D*f_bar(1:end/2),PF20(1:end/2),tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
            PF2_1 = PF2_1 + dot(S.W,f_bar(1:end/2))/sum(S.W);
            PF2_2=AAR(B,D*f_bar(end/2+1:end),PF20(end/2+1:end),tol,maxit,[],[],[],[],[]);
            PF2_2 = PF2_2 + dot(S.W,f_bar(end/2+1:end))/sum(S.W);
            PF2(:,k) = vertcat(PF2_1,PF2_2);
        end
        PF20 = PF2(:,k);
        if rem(iter,p)==0                          % apply every p iters        

            
            % Yk = (DF'*DF)\(DF'*f1);

            dp = DX*Yk + mix*PF20;
        else
            dp=0;
        end    
    end
    x_new = x1 + mix*PF0 ;%- dp;
    Xold = x1;
    Fold = f1;
else
    dp1 = 0;
    dp2 = 0;
    if iter == 1
        PF1=AAR(B,D*f1(:,1),[],tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
        PF0_1=PF1;
        PF1=AAR(B,D*f1(:,2),[],tol,maxit,[],[],[],[],[]);                % The 1st preconditioned term
        PF0_2=PF1;
    else
        PF1=AAR(B,D*f1(:,1),PF0_1,tol,maxit,[],[],[],[],[]); 
        PF0_1=PF1;
        PF1=AAR(B,D*f1(:,2),PF0_2,tol,maxit,[],[],[],[],[]); 
        PF0_2=PF1;
    end
    

    if iter > 1
        k = mod(iter-2,m)+1;

        DX1(:,k) = x1(:,1) - Xold(:,1);
        DF1(:,k) = f1(:,1) - Fold(:,1);
        DX2(:,k) = x1(:,2) - Xold(:,2);
        DF2(:,k) = f1(:,2) - Fold(:,2);
        
        Yk1 = pinv(DF1'*DF1)*(DF1'*f1(:,1));
        Yk2 = pinv(DF2'*DF2)*(DF2'*f1(:,2));
        f_bar1 = DF1*Yk1;
        f_bar2 = DF2*Yk2;
        if iter == 2
            PF2=AAR(B,D*f_bar1,[],tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
            PF20_1 = PF2;
            PF2=AAR(B,D*f_bar2,[],tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
            PF20_2 = PF2;
        else
            PF2=AAR(B,D*f_bar1,PF20_1,tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
            PF20_1 = PF2;
            PF2=AAR(B,D*f_bar2,PF20_2,tol,maxit,[],[],[],[],[]);  % The 2nd precondtioned term
            PF20_2 = PF2;
        end
        
        
        if rem(iter,p)==0                          % apply every p iters 
            dp1 = DX1*Yk1 + mix*PF20_1;
            dp2 = DX2*Yk2 + mix*PF20_2;
        else
            dp1 = 0;
            dp2 = 0;
        end    
    end
    x_new(:,1) = x1(:,1) + mix*PF0_1 ;%- dp1;
    x_new(:,2) = x1(:,2) + mix*PF0_2 ;%- dp2;

    Xold = x1;
    Fold = f1;
end
    
