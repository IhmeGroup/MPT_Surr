%% CREATE HYPERELLIPSOIDS/CUBOIDS FOR ALL CONSTRAINTS
% AUTHOR : G. PAVAN BHARADWAJ
% REFER TO LICENSE.pdf ON REPOSITORY FOR USAGE RESTRICTIONS

% clear all;
clc; clearvars;close all; yalmip('clear');cvx_clear;

%% INPUTS

% SURROGATE NAME
surr = 'violi';

% SHAPE (cuboid, ellipsoid, both)
shape = 'both';

save_output = 0;
no_graph = 1;
draw_2d = 1;
compute_distillation = 0;
linear_dcn = 1;

% VERBOSITY (1 for output)
% 2 for verification
verbose = 1;

% TARGET PROPERTIES FLAG
useCompositionOnly = true;

% RELATIVE ERROR THRESHOLDS
eps_M = 0.05;
eps_HC = 0.005;
eps_DC = 0.0125;
eps_IDT = 0.02;
eps_DCN = 0.02;
eps_TSI = 0.04;

% CURRENT CODE PATH
curr_code_path = pwd;

% MECHANISM PATH
mech_path = '/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Surrogates/Mechanisms/POLIMI_TOT';

% DISTILLATION CODE PATH
dist_code_path = '/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Multicomponent/Droplet_Evap/src';

% IDT DUMP PATH
idt_data_path = 'Data/won_eps/';

% DISTILLATION DATA PATH
dist_data_path = 'Data/won_eps/dc_err/';

%% LOAD MECHANISM
adddir(mech_path);
addpath('./include/');

if (~exist('gas'))
  gas = Solution('POLIMI_TOT.xml','gas');
end

disp('Loaded mechanism');

%% SET TARGET DETAILS
[palette,palette_label,palette_tsi,palette_dcn,exp_comp,target_mw,target_hc] = ...
    set_target_details(surr);

%% SET GLOBAL QUANTITIES

% SPECIES INDEX
palette_size = length(palette);
sp_index = zeros(1,palette_size);
for i = 1:palette_size
    sp_index(i) = speciesIndex(gas,palette{i});
end

assert(all(sp_index~=0));
disp('Identified palette compounds');

%% LOAD NONLINEAR DATA FROM pyIDT

file_name = strcat(idt_data_path,surr,'.dat');
data = dlmread(file_name,',');
%% COMPUTE DISTILLATION FILTER

if (compute_distillation)
  % SET DISTILLATION CODE PATH
  chdir(dist_code_path);
    
  err = zeros(size(data,1),1);
  disp(strcat('Number of points : ',num2str(size(data,1))));
  
  % COMPUTE DISTILLATION CURVE USING SEPARATE CODE
  for i = 1:size(data,1)
      switch(surr)
          case('dooley')
            compounds_file = 'dryer_surr';
          otherwise
            compounds_file = 'vihans_surr';
      end
      err(i) = distMoleFrac(compounds_file,false,data(i,1:palette_size),'posf4658_simdist','fractional');
      disp([i,data(i,1:palette_size),err(i)]);
  end
  
  % WRITE DC ERROR TO FILE
  chdir(curr_code_path);
  save(strcat(dist_data_path,surr,'_err.mat'),'err');

else

  % OTHERWISE, JUST LOAD DISTILLATION DATA
  file_name = strcat(dist_data_path,surr,'_err.mat');
  load(file_name,'err');

end

% FILTER USING DISTILLATION CURVE
[data, target_dc] = filter_dc(data,exp_comp,eps_DC,err);

% FILTER USING MW AND HC
data = filter_mw(data,exp_comp,eps_M,gas,sp_index);
data = filter_hc(data,exp_comp,eps_HC,gas,sp_index);

% FILTER USING TSI
[data,target_tsi] = filter_tsi(data,exp_comp,eps_TSI,palette_tsi);

% FILTER USING IDT/BLENDED DCN
if (~linear_dcn)
   [data, target_idt] = filter_idt(data,exp_comp,eps_IDT);
else
   [data, target_dcn] = filter_dcn(data,exp_comp,eps_DCN,palette_dcn);
end
%% CREATE CONVEX HULL

P = Polyhedron(data(:,1:palette_size-1));

% CHECK IF POLYHEDRON IS FULL DIMENSIONAL
if (P.isFullDim())
  disp('Polyhedron is full dimensional');
else
  disp('Polyhedron is NOT full dimensional');
end
    
%% ELLIPSOIDS

if (strcmp(shape,'ellipsoid')) 
  [Ell_AA,E_AA,Ell,E,x0_AA,x0,y0] = get_ellipsoids(P);

  if (~no_graph)
      if (draw_2d)
          draw_ellipsoids_2d(P,Ell,Ell_AA,exp_comp,palette_label,save_output);
      else
          draw_ellipsoids_3d(P,Ell,Ell_AA,exp_comp,palette_label);
      end
  end

  % PRINT STATISTICS
  if (verbose==1)
    disp(palette);
    disp('Exp. composition: ');
    disp(num2str(exp_comp));
    
    disp_ellipse_stats(P,x0_AA,E_AA,x0,E);
  end

elseif (strcmp(shape,'cuboid'))


  %% FIT HYPERCUBOIDS
  A = P.A;
  b = P.b;

 [x0_R,diff_R,y0] = getInnerCuboid(A,b);
 [x0_out,diff_out] = getOuterCuboid(A,b);

 %% PLOT HYPERCUBOIDS

 if (~no_graph)
     if (draw_2d)
         draw_cuboids_2d(P,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output);
     else
         draw_cuboids_3d(P,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label);
     end
 end
  
  % PRINT STATISTICS
  if (verbose==1)
    disp(palette);
    disp('Exp. composition: ');
    disp(num2str(exp_comp));
    
    disp_cuboid_stats(P,x0_R,diff_R,x0_out,diff_out,y0);
  end

elseif (strcmp(shape,'both'))
    
  % GET ELLIPSOIDS  
  [Ell_AA,E_AA,Ell,E,x0_AA,x0,y0] = get_ellipsoids(P);
  
  % FIT HYPERCUBOIDS
  A = P.A;
  b = P.b;

 [x0_R,diff_R,y0] = getInnerCuboid(A,b);
 [x0_out,diff_out] = getOuterCuboid(A,b);
 
 if (~no_graph)
     if (draw_2d)
         draw_both_2d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output,surr);
     else
         draw_both_3d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label);
     end
 end
  
 % PRINT STATISTICS
  if (verbose==1)
    disp(palette);
    disp('Exp. composition: ');
    disp(num2str(exp_comp));
    
    disp_ellipse_stats(P,x0_AA,E_AA,x0,E);
    disp_cuboid_stats(P,x0_R,diff_R,x0_out,diff_out,y0);
    
  elseif (verbose == 2)
      
     % RELOAD NONLINEAR DATA
     chdir(curr_code_path);
     file_name = strcat(idt_data_path,surr,'.dat');
     data = dlmread(file_name,',');
     
     % ISOLATE IDT
     idts = abs(data(:,palette_size+1));
     
     % MW DETAILS
     MW_Vec_full = molecularWeights(gas);
     MW_Vec = zeros(palette_size,1);
     for i = 1:palette_size
         MW_Vec(i) = MW_Vec_full(sp_index(i));
     end
     
     % HC DETAILS
     H_Vec = zeros(palette_size,1);
     C_Vec = zeros(palette_size,1);
     for i = 1:palette_size
         H_Vec(i) = nAtoms(gas,sp_index(i),'H');
         C_Vec(i) = nAtoms(gas,sp_index(i),'C');
     end

     % Constraint bounds
     fprintf('MW : %f %f \n',target_mw*(1-eps_M),target_mw*(1+eps_M));
     fprintf('HC : %f %f \n',target_hc*(1-eps_HC),target_hc*(1+eps_HC));
     fprintf('TSI : %f %f \n',target_tsi*(1-eps_TSI),target_tsi*(1+eps_TSI));
     fprintf('IDT : %f %f \n',(10^-(target_idt*(1+eps_IDT)))),(10^-(target_idt*(1-eps_IDT)));
     fprintf('DC : %f %f \n',target_dc*(1-eps_DC),target_dc*(1+eps_DC));
     
     % INNER HYPERCUBOID ERRORS
     
     fprintf('\n INNER HYPERCUBOID \n');
     x_test = [x0_R,1-sum(x0_R)];
     
     % CALCULATE NEAREST COMPOSITION
     dist =  zeros(size(data,1),1);
     for i = 1:size(data,1)
         dist(i) = norm(x_test(1:palette_size-1)-data(i,1:palette_size-1));
     end
     [min_norm,nearest_index] = min(dist);
     
     % Compute
     fprintf('MW : %f \n', dot(x_test,MW_Vec));
     fprintf('HC : %f \n',dot(x_test,H_Vec)/dot(x_test,C_Vec));
     fprintf('TSI : %f \n',dot(x_test,palette_tsi));
     
     fprintf('IDT : %f \n',idts(nearest_index));
     fprintf('DC : %f \n',err(nearest_index));

     % AXIS-ALIGNED HYPERELLIPSOID ERRORS
 
     fprintf('\n AXIS-ALIGNED HYPERELLIPSOID \n');
     x_test = [x0_AA,1-sum(x0_AA)];
     
     % CALCULATE NEAREST COMPOSITION
     dist =  zeros(size(data,1),1);
     for i = 1:size(data,1)
         dist(i) = norm(x_test(1:palette_size-1)-data(i,1:palette_size-1));
     end
     [min_norm,nearest_index] = min(dist);
     
     % Compute
     fprintf('MW : %f \n', dot(x_test,MW_Vec));
     fprintf('HC : %f \n',dot(x_test,H_Vec)/dot(x_test,C_Vec));
     fprintf('TSI : %f \n',dot(x_test,palette_tsi));
     
     fprintf('IDT : %f \n',idts(nearest_index));
     fprintf('DC : %f \n',err(nearest_index));
     
     % MAXIMUM VOLUME HYPERELLIPSOID ERRORS
     
     fprintf('\n MAXIMUM-VOLUME HYPERELLIPSOID \n');
     x_test = [x0,1-sum(x0)];
     
     % CALCULATE NEAREST COMPOSITION
     dist =  zeros(size(data,1),1);
     for i = 1:size(data,1)
         dist(i) = norm(x_test(1:palette_size-1)-data(i,1:palette_size-1));
     end
     [min_norm,nearest_index] = min(dist);
     
     % Compute
     fprintf('MW : %f \n', dot(x_test,MW_Vec));
     fprintf('HC : %f \n',dot(x_test,H_Vec)/dot(x_test,C_Vec));
     fprintf('TSI : %f \n',dot(x_test,palette_tsi));
     
     fprintf('IDT : %f \n',idts(nearest_index));
     fprintf('DC : %f \n',err(nearest_index));
       
  end
 
  if (save_output || no_graph)
      close all;
  end

end