function [data_out,target_tsi] = filter_tsi(data,exp_comp,eps,palette_tsi)

palette_size = size(data,2)-1;

% SET TARGET PROPERTIES
target_tsi = exp_comp*palette_tsi;

% REMOVE USING TSI
lb = (1 - eps)*target_tsi;
ub = (1 + eps)*target_tsi;

% OBTAIN MW
tsi =  zeros(size(data,1),1);
for i = 1:size(data,1)
    tsi(i) = dot(data(i,1:palette_size),palette_tsi);
end

data_out = data(tsi <= ub & tsi >= lb,:);

disp('Removed using TSI');
disp(strcat('Number of points: ',num2str(size(data_out,1))));

end