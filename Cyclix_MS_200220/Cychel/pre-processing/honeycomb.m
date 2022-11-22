a = 2.67; % in bohr
A = sqrt(3)*a;
a1 = [sqrt(3)/2 1/2]*A; a2 = [sqrt(3)/2 -1/2]*A;
D1 = 0; D2 = 0.5;
x_basis = [0 0 D1; horzcat(a1/3+a2/3,D2)];

natoms1 = 225;
natoms2 = 240;
x_lattice1 = zeros(natoms1,3);
x_lattice2 = zeros(natoms2,3);
count1 = 1;
count2 = 1;
for i = 0:14
	for j = -i:i
		x_lattice1(count1,:) = x_basis(1,:) + horzcat(i*a1 + j*a2,D1);
		x_lattice2(count2,:) = x_basis(2,:) + horzcat(i*a1 + j*a2,D2);
		count1 = count1+1;
		count2 = count2+1;
	end
	x_lattice2(count2,:) = x_basis(2,:) + horzcat(i*a1 + (-i-1)*a2,D2);
	count2 = count2+1;
end

x_lattice1 = x_lattice1(x_lattice1(:,1) < 62,:,:);
x_lattice1 = x_lattice1(x_lattice1(:,2) < 33,:,:);

x_lattice2 = x_lattice2(x_lattice2(:,1) < 62,:,:);
x_lattice2 = x_lattice2(x_lattice2(:,2) < 33,:,:);

figure()
hold on
scatter(x_lattice1(:,1),x_lattice1(:,2));
scatter(x_lattice2(:,1),x_lattice2(:,2));
hold off
saveas(gcf,'honeycomb','epsc')
%file = strcat('C','.mat');
%file_loc = fullfile('./',file);
%save(file_loc,'x_lattice');

x_lattice1(:,2) = 0;
x_lattice2(:,2) = 0;
fid = fopen('Graphene.xyz', 'w');
fprintf(fid, '%d\n',2*size(x_lattice1,1));
fprintf(fid, 'Graphene\n');
for i = 1:size(x_lattice1,1)
	fprintf(fid,'C %f %f %f\n',x_lattice1(i,:));
	fprintf(fid,'O %f %f %f\n',x_lattice2(i,:));
end
fclose(fid)
