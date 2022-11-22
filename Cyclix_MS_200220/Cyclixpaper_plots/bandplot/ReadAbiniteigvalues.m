function ReadAbiniteigvalues(element)

filename = strcat('./Abinit_eigs/',element,'_EIG');
[fid,~] = fopen(filename,'r');  
assert(fid~=-1,'Error: Check the element name and the corresponding _EIG file') ;
fscanf(fid,'%s',4); 
nkpt = fscanf(fid,'%d',1);
textscan(fid,'%s',1,'delimiter','\n');
fscanf(fid,'%s',3);
nband = fscanf(fid,'%d',1);
E_abinit = zeros(nkpt,nband);

for kpt = 1:nkpt
    textscan(fid,'%s',1,'delimiter','\n');
    E_abinit(kpt,:) = fscanf(fid,'%g');% expected to be in Ha
end

file = strcat('ABINITEigenvalues_',element,'.mat');
save(file,'E_abinit') ;
    
    
