function [data_out,target_dcn] = filter_dcn(data,exp_comp,eps_DCN,palette_dcn)

palette_size = size(data,2)-1;

% SET TARGET PROPERTIES
target_dcn = exp_comp*palette_dcn;

% REMOVE USING MW
lb = (1 - eps_DCN)*target_dcn;
ub = (1 + eps_DCN)*target_dcn;

% OBTAIN DCN
dcn =  zeros(size(data,1),1);
for i = 1:size(data,1)
    dcn(i) = dot(data(i,1:palette_size),palette_dcn);
end

data_out = data(dcn <= ub & dcn >= lb,:);

disp('Removed using DCN');
disp(strcat('Number of points: ',num2str(size(data_out,1))));

end