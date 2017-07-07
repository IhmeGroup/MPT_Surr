function [data_out,target] = filter_dc(data,exp_comp,eps,err)

palette_size = size(data,2)-1;

% FIND CLOSEST COMPOSITION TO EXPERIMENT
dist =  zeros(size(data,1),1);
for i = 1:size(data,1)
    dist(i) = norm(exp_comp(1:palette_size-1)-data(i,1:palette_size-1));
end
[min_norm,nearest_index] = min(dist);

target = err(nearest_index);
lb = (1 - eps)*target;
ub = (1 + eps)*target;

data_out = data(err <= ub & err >= lb,:);

disp('Removed using DC');
disp(strcat('Number of points left: ',num2str(size(data_out,1))));

end