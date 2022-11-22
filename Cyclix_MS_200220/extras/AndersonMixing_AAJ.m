function [x_new] = AndersonMixing_AAJ(x1, f1, m, iter, mix, p, spin)
% f = residual, m = history, x=iterate, f(x)=g(x)-x, p=AJ step
% Reference:
%         Pratapa, Phanisri P., Phanish Suryanarayana, and John E. Pask. 
%         "Anderson acceleration of the Jacobi iterative method: 
%         An efficient alternative to Krylov methods for large, sparse linear systems." 
%         Journal of Computational Physics 306 (2016): 43-54.
persistent DX1 DX2 DX DF1 DF2 DF Xold Fold 
if spin == 1
    dp = 0;
    if iter > 1
        k = mod(iter-2,m)+1;

        DX(:,k) = x1 - Xold;
        DF(:,k) = f1 - Fold;

        if rem(iter,p)==0        % apply every p iters        

            Yk = pinv(DF'*DF)*(DF'*f1);        
            % Yk = (DF'*DF)\(DF'*f1);

            dp = (DX + mix*DF)*Yk;
        else
            dp=0;
        end    
    end
    x_new = x1 + mix*f1 - dp;

    Xold = x1;
    Fold = f1;
elseif spin == 2
    
    %dp = 0;
    dp1 = 0;
    dp2 = 0;
    if iter > 1
        k = mod(iter-2,m)+1;

        DX1(:,k) = x1(:,2) - Xold(:,2);
        DF1(:,k) = f1(:,2) - Fold(:,2);
        
        DX2(:,k) = x1(:,3) - Xold(:,3);
        DF2(:,k) = f1(:,3) - Fold(:,3);
        
        %DX(:,k) = x1(:,1) - Xold(:,1);
        DF(:,k) = f1(:,1) - Fold(:,1);

        if rem(iter,p)==0        % apply every p iters        

            Yk = pinv(DF'*DF)*(DF'*f1(:,1));       
            % Yk = (DF'*DF)\(DF'*f1);

            %dp = (DX + mix*DF)*Yk;
            dp1 = (DX1 + mix*DF1)*Yk;
            dp2 = (DX2 + mix*DF2)*Yk;


%             dp1 = (DX1 + mix*DF)*Yk;
%             dp2 = (DX2 + mix*DF)*Yk;
        else
            %dp=0;
            dp1=0;
            dp2=0;
        end    
    end
    x_new(:,2) = x1(:,2) + mix*f1(:,2) - dp1;
    x_new(:,3) = x1(:,3) + mix*f1(:,3) - dp2;
    x_new(:,1) = x_new(:,2) + x_new(:,3);

%     x_new(:,2) = x1(:,2) + mix*f1(:,1) - dp;
%     x_new(:,3) = x1(:,3) + mix*f1(:,1) - dp;
%     x_new(:,1) = x_new(:,2) + x_new(:,3);

%     x_new(:,2) = x1(:,2) + mix*f1(:,1) - dp1;
%     x_new(:,3) = x1(:,3) + mix*f1(:,1) - dp2;
%     x_new(:,1) = x_new(:,2) + x_new(:,3);
    

    Xold = x1;
    Fold = f1;
elseif spin == 3
    dp1 = 0;
    dp2 = 0;
    if iter > 1
        k = mod(iter-2,m)+1;

        DX1(:,k) = x1(:,1) - Xold(:,1);
        DF1(:,k) = f1(:,1) - Fold(:,1);
        
        DX2(:,k) = x1(:,2) - Xold(:,2);
        DF2(:,k) = f1(:,2) - Fold(:,2);
        
        if rem(iter,p)==0        % apply every p iters        

            Yk1 = pinv(DF1'*DF1)*(DF1'*f1(:,1));
            Yk2 = pinv(DF2'*DF2)*(DF2'*f1(:,2));
            
            dp1 = (DX1 + mix*DF1)*Yk1;
            dp2 = (DX2 + mix*DF2)*Yk2;
        else
            dp1=0;
            dp2=0;
        end    
    end
    x_new(:,1) = x1(:,1) + mix*f1(:,1) - dp1;
    x_new(:,2) = x1(:,2) + mix*f1(:,2) - dp2;

    Xold = x1;
    Fold = f1;
end

