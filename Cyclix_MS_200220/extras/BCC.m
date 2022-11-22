a = [0 0 0;0.5 0.5 0.5];
ncells = [5 5 5];
ax = zeros(size(a,1)*ncells(1),3);
b = a;
b(:,1) = b(:,1)/ncells(1);
for i = 1:ncells(1)
    ax(2*i-1:2*i,:) = b + [(i-1)/ncells(1) 0 0;(i-1)/ncells(1) 0 0];
end

ay = repmat(ax,ncells(2),1);
step = size(ax,1);
for i = 1:ncells(2)  
    ay(step*(i-1)+1:step*i,2) = ax(:,2)/ncells(2) + (i-1)/ncells(2);
end

az = repmat(ay,ncells(3),1);
step = size(ay,1);
for i = 1:ncells(3)  
    az(step*(i-1)+1:step*i,3) = ay(:,3)/ncells(3) + (i-1)/ncells(3);
end