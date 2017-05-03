function [MW_Vec,H_Vec,C_Vec,A,b] = linear_constraints(gas,exp_comp,useCompositionOnly)

global palette_size sp_index;

% MW DETAILS
MW_Vec_full = molecularWeights(gas);
MW_Vec = zeros(palette_size,1);
for i = 1:palette_size
    MW_Vec(i) = MW_Vec_full(sp_index(i));
end

MW_Vec_T = MW_Vec(1:end-1);
MW_Vec_L = MW_Vec(end);

% HC DETAILS
H_Vec = zeros(palette_size,1);
C_Vec = zeros(palette_size,1);
for i = 1:palette_size
   H_Vec(i) = nAtoms(gas,sp_index(i),'H');
   C_Vec(i) = nAtoms(gas,sp_index(i),'C');
end

H_Vec_T = H_Vec(1:end-1);
H_Vec_L = H_Vec(end);

C_Vec_T = C_Vec(1:end-1);
C_Vec_L = C_Vec(end);

% SET TARGET PROPERTIES
if (useCompositionOnly)
  target_mw = exp_comp*MW_Vec;
  target_hc = (exp_comp*H_Vec)/(exp_comp*C_Vec);
end

% CONSTRUCT MW MATRIX
b_MW = [target_mw*(1 + eps_M);-target_mw*(1 - eps_M)];
A_MW = [MW_Vec';-MW_Vec'];

% CONSTRUCT HC MATRIX
b_HC = [0;0];
A_HC = [H_Vec' - C_Vec'*(target_hc*(1 + eps_HC)); ...
       -(H_Vec' - C_Vec'*(target_hc*(1 - eps_HC)))];

% CONSTRUCT TRUNCATED MW MATRIX
b_MW_T = [target_mw*(1 + eps_M)-MW_Vec_L; -target_mw*(1 - eps_M)+MW_Vec_L];
A_MW_T = [MW_Vec_T'-MW_Vec_L;-MW_Vec_T'+MW_Vec_L];

% CONSTRUCT TRUNCATED HC MATRIX
b_HC_T = [- H_Vec_L + C_Vec_L*(target_hc*(1 + eps_HC));...
        -(-H_Vec_L + C_Vec_L*(target_hc*(1 - eps_HC)))];
A_HC_T = [(H_Vec_T'- H_Vec_L) - (C_Vec_T'-C_Vec_L)*(target_hc*(1 + eps_HC)); ...
       -((H_Vec_T'-H_Vec_L) - (C_Vec_T'-C_Vec_L)*(target_hc*(1 - eps_HC)))];
   
A = [A_MW_T;A_HC_T];
b = [b_MW_T;b_HC_T];

% CREATE STANDARD CONSTRAINTS

% POSITIVITY AND LEQ 1
A = [A;-eye(palette_size-1);eye(palette_size-1)];
b = [b;zeros(palette_size-1,1);ones(palette_size-1,1)];

end