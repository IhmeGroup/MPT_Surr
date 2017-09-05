%% CREATE HYPERELLIPSOIDS/CUBOIDS FOR MW AND HC RATIO CONSTRAINTS
% AUTHOR : G. PAVAN BHARADWAJ
% REFER TO LICENSE.pdf ON REPOSITORY FOR USAGE RESTRICTIONS

% clear all;
clc; clearvars; close all; yalmip('clear');cvx_clear;

%% INPUTS

% SURROGATE NAME
surr = 'dooley';
shape = 'cuboid';
save_output = 0;
no_graph = 0;

% DRAW 2D OR 3D PROJECTIONS (0 for 3D and 1 for 2D)
draw_2d = 1;

% VERBOSITY (1 for output)
verbose = 1;

% TARGET PROPERTIES FLAG
useCompositionOnly = true;

% RELATIVE ERROR THRESHOLDS
eps_M = 0.05;
eps_HC = 0.005;

% MECHANISM PATH
mech_path = '/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Surrogates/Mechanisms/POLIMI_TOT';

%% LOAD MECHANISM
adddir(mech_path);
addpath('./include/');

if (~exist('gas'))
  gas = Solution('POLIMI_TOT.xml','gas');
end

disp('Loaded mechanism');

%% SET TARGET DETAILS
[palette,palette_label,~,~,exp_comp,target_mw,target_hc] = ...
    set_target_details(surr);

%%  SET GLOBAL QUANTITIES

% SPECIES INDEX
palette_size = length(palette);
sp_index = zeros(1,palette_size);
for i = 1:palette_size
    sp_index(i) = speciesIndex(gas,palette{i});
end

assert(all(sp_index~=0));
disp('Identified palette compounds');

%% CREATE LINEAR CONSTRAINT MATRIX

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

disp('Linear constraints');

%% CREATE POLYHEDRON

% CREATE POLYTOPE
P = Polyhedron(A,b);

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
 [x0_R,diff_R,y0] = getInnerCuboid(A,b);
 disp('Solving outer cuboid program...');
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
  
elseif(strcmp(shape,'both'))
  % GET ELLIPSOIDS  
  [Ell_AA,E_AA,Ell,E,x0_AA,x0,y0] = get_ellipsoids(P);
  
  % FIT HYPERCUBOIDS
  A = P.A;
  b = P.b;

 [x0_R,diff_R,y0] = getInnerCuboid(A,b);
 [x0_out,diff_out] = getOuterCuboid(A,b);
 
 if (~no_graph)
     if (draw_2d)
         draw_both_2d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output);
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
  end
 
  if (save_output || no_graph)
      close all;
  end
    
end

if (save_output)
  close all;
end
