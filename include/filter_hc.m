function [data_out] = filter_hc(data,exp_comp,eps,gas,sp_index)

palette_size = size(data,2)-1;

% HC DETAILS
H_Vec = zeros(palette_size,1);
C_Vec = zeros(palette_size,1);
for i = 1:palette_size
   H_Vec(i) = nAtoms(gas,sp_index(i),'H');
   C_Vec(i) = nAtoms(gas,sp_index(i),'C');
end

% SET TARGET PROPERTIES
target_hc = (exp_comp*H_Vec)/(exp_comp*C_Vec);

% REMOVE USING HC
lb = (1 - eps)*target_hc;
ub = (1 + eps)*target_hc;

% OBTAIN HC RATIO
hc =  zeros(size(data,1),1);
for i = 1:size(data,1)
    vec = data(i,1:palette_size);
    hc(i) = dot(vec,H_Vec)/dot(vec,C_Vec);
end

data_out = data(hc <= ub & hc >= lb,:);

disp('Removed using HC');
disp(strcat('Number of points left: ',num2str(size(data_out,1))));

end