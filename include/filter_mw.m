function [data_out] = filter_mw(data,exp_comp,eps_M,gas,sp_index)

palette_size = size(data,2)-1;

% MW DETAILS
MW_Vec_full = molecularWeights(gas);
MW_Vec = zeros(palette_size,1);
for i = 1:palette_size
    MW_Vec(i) = MW_Vec_full(sp_index(i));
end

% SET TARGET PROPERTIES
target_mw = exp_comp*MW_Vec;

% REMOVE USING MW
lb = (1 - eps_M)*target_mw;
ub = (1 + eps_M)*target_mw;

% OBTAIN MW
mw =  zeros(size(data,1),1);
for i = 1:size(data,1)
    mw(i) = dot(data(i,1:palette_size),MW_Vec);
end

data_out = data(mw <= ub & mw >= lb,:);

disp('Removed using MW');
disp(strcat('Number of points: ',num2str(size(data_out,1))));

end