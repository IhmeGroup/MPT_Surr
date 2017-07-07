function [data_out,target] = filter_idt(data,exp_comp,eps)

palette_size = size(data,2)-1;

% ISOLATE IDT
idts = abs(log10(data(:,palette_size+1)));

% FIND CLOSEST COMPOSITION TO EXPERIMENT
dist =  zeros(size(data,1),1);
for i = 1:size(data,1)
    dist(i) = norm(exp_comp(1:palette_size-1)-data(i,1:palette_size-1));
end
[min_norm,nearest_index] = min(dist);

target = idts(nearest_index);
lb = (1 - eps)*target;
ub = (1 + eps)*target;

data_out = data(idts <= ub & idts >= lb,:);

disp('Removed using IDT');
disp(strcat('Number of points left: ',num2str(size(data_out,1))));


end