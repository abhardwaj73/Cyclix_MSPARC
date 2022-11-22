% GENERATES DIAMOND CUBIC

format long;
clear all
clc

%%%%%%%%%%%%%% INPUT PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=1;   % unit cells in x direction
Ny=1;   % unit cells in y direction
Nz=1;   % unit cells in z direction
a = 7.212489168102784;  % lattice constant (length of unit cell)
aa=a/2;

Rx = a * Nx;
Ry = a * Ny;
Rz = a * Nz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%          Basis Vectors of lattice                %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%  Diamond cubic   %%%%%%%%%
b1 = aa*[1 1 0] ;
b2 = aa*[1 0 1] ;
b3 = aa*[0 1 1] ;
shift = aa/2*[1 1 1] ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%              Supercell information                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mm = 10 ;
count = 1 ;
%xa_sc = zeros(18,3) ;
tol = 1e-8 ;
for p1=-mm:mm
    for p2=-mm:mm
        for p3=-mm:mm
            temp = p1*b1 + p2*b2 + p3*b3 ; 
            d = norm(temp);
            if (temp(1)>=-(tol) && temp(1)<=(Nx*a+tol) ...
                    && temp(2)>=-(tol) && temp(2)<=(Ny*a+tol) ...
                    && temp(3)>=-(tol) && temp(3)<=(Nz*a+tol))
                %store positions of accepted atoms
                pos_Si(count,:) = temp; 
                %Store the basis vectors of accepted atoms                
                A(count,1:3) = [p1,p2,p3];                
                count = count + 1;                
            end
            % temp = p1*b1 + p2*b2 + p3*b3 + shift ;
            temp = temp + shift ;
            d = norm(temp);
            if (temp(1)>=-(tol) && temp(1)<=(Nx*a+tol) ...
                    && temp(2)>=-(tol) && temp(2)<=(Ny*a+tol) ...
                    && temp(3)>=-(tol) && temp(3)<=(Nz*a+tol))
                %store positions of accepted atoms
                pos_Si(count,:) = temp;
                %Store the basis vectors of accepted atoms                
                A(count,1:3) = [p1,p2,p3];                
                count = count + 1;           
            end
        end
    end
end
disp('Number of Silicon Atoms :');
nSi = size(pos_Si,1);
disp(nSi);
%plot3(pos_Si(:,1),pos_Si(:,2),pos_Si(:,3),'ro')

% need to remove face atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Removing atoms from 3 faces due to PBC's %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol1 = 1e-6;
ind1 = ones(size(pos_Si,1),1);
for qq=1:size(pos_Si,1)    
    if abs(pos_Si(qq,1)-Rx)<tol1 || abs(pos_Si(qq,2)-Ry)<tol1 || abs(pos_Si(qq,3)-Rz)<tol1
       ind1(qq) = 0;
    end    
end
ind1 = ind1 > 0.5 ;
pos_Si = pos_Si(ind1,:);



%% plot atom positions
plot3(pos_Si(:,1),pos_Si(:,2),pos_Si(:,3),'o','markersize',20, ...
    'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0],'linewidth',1);

%% mark atom coordinates
for k=1:size(pos_Si,1)
      text(pos_Si(k,1)+0.25,pos_Si(k,2)+0.25,pos_Si(k,3)+0.25,['(' num2str(pos_Si(k,1)) ',' ...
          num2str(pos_Si(k,2)) ',' num2str(pos_Si(k,3)) ')'],'fontsize',16)
end
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
set(gca,'FontSize',16)
% axis equal

%% draw cuboids
% vert = [0 0 0;1 0 0;1 1 0;0 1 0;0 0 1;1 0 1;1 1 1;0 1 1];
% fac = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8];
% for k = 1:Nx*Ny*Nz 
    % kk_e = ceil(k/(Nx*Ny));
    % jj_e = ceil( (k-(kk_e-1)*Nx*Ny)/Nx );
    % ii_e = k - (kk_e-1)*Nx*Ny - (jj_e-1)*Nx;
    % fprintf('k = %d, %d %d %d\n',k,ii_e,jj_e,kk_e);
    % vert_k = vert + kron(ones(8,1),[ii_e jj_e kk_e]-1);
    % vert_k = a * vert_k;
    % patch('Vertices',vert_k,'Faces',fac,'FaceVertexCData',hsv(8),'FaceColor','none','linewidth',2)
% end
    
% %% view set up    
% axis equal
% view([-1 -1.5 0.7])

%% print atom positions
pos_Si
