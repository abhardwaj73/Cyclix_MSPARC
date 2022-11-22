d = 13;
twist = 0.072894889164856;
tau = 4.0332045685;
ncell_c = d;
ncell_z = 11;

xcart = [8.617920886117433    4.133073094220444    0.676447793032692
         8.445317000224508    4.475915682699788    3.341240136113001
];

natom_u = size(xcart,1);
xcyc = zeros(natom_u*ncell_c*ncell_z,3);
xcyc(1:natom_u,:) = [sqrt(xcart(:,1).^2+xcart(:,2).^2)  atan2(xcart(:,2),xcart(:,1))  xcart(:,3)];

gen_c = [0 2*pi/d 0];
Psi = twist*tau;
gen_z = [0 Psi tau];

for i = 1:ncell_c-1
    xcyc(natom_u*i+1:natom_u*(i+1),:) = xcyc(1:natom_u,:) + gen_c*i;
end

for i = 1:ncell_z-1
    xcyc(natom_u*ncell_c*i+1:natom_u*ncell_c*(i+1),:) = xcyc(1:natom_u*ncell_c,:) + gen_z*i;
end

fprintf('\ncylindrical coordinates\n');
xcyc

xcart_full = [xcyc(:,1).*cos(xcyc(:,2)) xcyc(:,1).*sin(xcyc(:,2)) xcyc(:,3)]

fid = fopen('Nanotube.xyz', 'w');
fprintf(fid, '%d\n',size(xcart_full,1));
fprintf(fid, 'Nanotube\n');
for i = 1:2:size(xcart_full,1)
    fprintf(fid,'C %f %f %f\n',xcart_full(i,:));
    fprintf(fid,'O %f %f %f\n',xcart_full(i+1,:));
end
fclose(fid);

