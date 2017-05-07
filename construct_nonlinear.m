%% CREATE HYPERELLIPSOIDS/CUBOIDS FOR ALL CONSTRAINTS
% AUTHOR : G. PAVAN BHARADWAJ
% REFER TO LICENSE.pdf ON REPOSITORY FOR USAGE RESTRICTIONS

% clear all;
clc; close all; yalmip('clear');cvx_clear;

%% INPUTS

% SURROGATE NAME
surr = 'hanson_a';

% SHAPE (cuboid, ellipsoid, both)
shape = 'both';

save_output = 0;
no_graph = 0;
draw_2d = 0;
compute_distillation = 0;

% VERBOSITY (1 for output)
verbose = 1;

% TARGET PROPERTIES FLAG
useCompositionOnly = true;

% RELATIVE ERROR THRESHOLDS
eps_M = 0.125;
eps_HC = 0.125;
eps_DC = 0.125;
eps_IDT = 0.075;
eps_TSI = 0.250;

% MECHANISM PATH
mech_path = '/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Surrogates/Mechanisms/POLIMI_TOT';

% DISTILLATION CODE PATH
dist_code_path = '/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Multicomponent/Droplet_Evap/src';

% IDT DUMP PATH
idt_data_path = 'Data/eps_half/CSV/';

% DISTILLATION DATA PATH
dist_data_path = 'Data/eps_half/dc_err/';

%% LOAD MECHANISM
adddir(mech_path);
addpath('./include/');

if (~exist('gas'))
  gas = Solution('POLIMI_TOT.xml','gas');
end

disp('Loaded mechanism');

%% SET TARGET DETAILS
[palette,palette_label,palette_tsi,exp_comp,target_mw,target_hc] = ...
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

file_name = strcat(idt_data_path,surr,'.out');
data = csvread(file_name);
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
      err(i) = distMoleFrac_obj(compounds_file,data(i,1:palette_size),'posf4658_simdist');
      disp([i,data(i,1:palette_size),err(i)]);
  end

else

  % OTHERWISE, JUST LOAD DISTILLATION DATA
  file_name = strcat(dist_data_path,surr,'_err.mat');
  load(file_name,'err');

end

% FILTER USING DISTILLATION CURVE
data = filter_dc(data,exp_comp,eps_DC,err);

% FILTER USING MW AND HC
data = filter_mw(data,exp_comp,eps_M,gas,sp_index);
data = filter_hc(data,exp_comp,eps_HC,gas,sp_index);

% FILTER USING TSI
data = filter_tsi(data,exp_comp,eps_TSI,palette_tsi);

% FILTER USING IDT
data = filter_idt(data,exp_comp,eps_IDT);
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

  if (draw_2d)
    draw_ellipsoids_2d(P,Ell,Ell_AA,exp_comp,palette_label,save_output);
  else
    draw_ellipsoids_3d(P,Ell,Ell_AA,exp_comp,palette_label);
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

  if (draw_2d)
    draw_cuboids_2d(P,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output);
  else
    draw_cuboids_3d(P,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label);
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
 
 if (draw_2d)
    draw_both_2d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output);
  else
    draw_both_3d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label);
  end
 
end

if (save_output || no_graph)
  close all;
end

