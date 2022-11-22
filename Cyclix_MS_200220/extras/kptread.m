function kptread(tnkpt)
fid1=fopen('k2','r'); 
kptgrid = zeros(tnkpt,3);
wkpt = zeros(tnkpt,1);
for i = 1:tnkpt
kptgrid(i,:) = cell2mat(textscan(fid1,'%f %f %f',1));
wkpt(i) = cell2mat(textscan(fid1,'%f',1,'delimiter','\n'));
end

save ('kptgrid','kptgrid');
save('wkpt','wkpt');
