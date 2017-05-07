%% SOLVE MPT PROGRAM FOR MW AND HC RATIO
% AUTHOR : G. PAVAN BHARADWAJ

clear all;
clc; close all; yalmip('clear');cvx_clear;

%% INPUTS

% SURROGATE NAME
surr = 'violi';

% TARGET PROPERTIES FLAG
useCompositionOnly = true;

%% LOAD MECHANISM
adddir('/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Surrogates/Mechanisms/POLIMI_TOT');
addpath('./include/');

if (~exist('gas'))
  gas = Solution('POLIMI_TOT.xml','gas');
end

disp('Loaded mechanism');

%% SET TARGET DETAILS
[palette,palette_label,~,exp_comp,target_mw,target_hc] = ...
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

% CONSTRUCT MW CONSTRAINT
MW_Vec_full = molecularWeights(gas);
MW_Vec = zeros(palette_size,1);
for i = 1:palette_size
    MW_Vec(i) = MW_Vec_full(sp_index(i));
end

% CONSTRUCT HC CONSTRAINT (NON-LINEAR)
H_Vec = zeros(palette_size,1);
C_Vec = zeros(palette_size,1);

for i = 1:palette_size
   H_Vec(i) = nAtoms(gas,sp_index(i),'H');
   C_Vec(i) = nAtoms(gas,sp_index(i),'C');
end

% SET TARGET PROPERTIES
if (useCompositionOnly)
  target_mw = exp_comp*MW_Vec;
  target_hc = (exp_comp*H_Vec)/(exp_comp*C_Vec);
end

%% MPT PROBLEM PARAMETERS

% ERROR THRESHOLDS
MW_err = 25.0;

% MPT VARIABLES 
x = sdpvar(palette_size, 1);
err = sdpvar(1, 1);

% CONSTRUCT CONSTRAINT MATRIX
b = [target_mw + err(1); -target_mw + err(1);0;0];
A = [MW_Vec';-MW_Vec';H_Vec' - target_hc*C_Vec'; -(H_Vec' - target_hc*C_Vec')];

% CONSTRAINTS
C = [A*x <= b, x >= 0, sum(x)==1, 0 <= err(1) <= MW_err];

%% SOLVE MPT PROBLEM
x_ax = linspace(0,MW_err,50);
data = [x_ax];

for i = 1:palette_size

    % OBJECTIVE FUNCTION
    sel = zeros(1,palette_size);
    sel(i) = 1;
    J = sel*x;
    
    % SOLVE MINIMUM PROBLEM
    plp = Opt(C, J, err, x);
    solution = plp.solve();
    
    for j = 1:length(x_ax)
        Z(j) = solution.xopt.feval(x_ax(j),'obj');
    end
    
    data = [data;Z];
    
    sel = zeros(1,palette_size);
    sel(i) = -1;
    J = sel*x;
    
    % SOLVE MAXIMUM PROBLEM
    plp = Opt(C, J, err, x);
    solution = plp.solve();
    
    for j = 1:length(x_ax)
        Z(j) = -solution.xopt.feval(x_ax(j),'obj');
    end
    
    data = [data;Z];

end

%% PLOT

figure('Position', [100, 100, 1000, 895]);
set(0,'DefaultAxesColorOrder',brewermap(palette_size,'Set1')); 

hold all;
% LOWER BOUNDS
for j = 1:palette_size
    plot(data(1,:),data(2*j,:),'--','LineWidth',2);
end

set(gca,'ColorOrderIndex',1);
% UPPER BOUNDS
k = zeros(1,palette_size);
for j = 1:palette_size
    k(j) = plot(data(1,:),data(2*j+1,:),'-','LineWidth',2);
end

set(gca,'ColorOrderIndex',1);
% EXPERIMENTAL COMPOSITION
for j = 1:palette_size
  h = plot(0,exp_comp(j),'^','MarkerSize',20);
  set(h,'MarkerEdgeColor',get(h,'Color'),'MarkerFaceColor',get(h,'Color'));
end

set(gcf,'Color',[1 1 1]);
set(gca,'FontName','Times New Roman','FontSize',32);
xlabel('\epsilon_{MW}','Interpreter','Tex');
ylabel('Mole Fraction');
h0 = plot(-1,-1,'^','MarkerSize',20,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
h1 = plot(0,0,'--k','LineWidth',2);
h2 = plot(0,0,'-k','LineWidth',2);
legend_array = [{'Experimental'};{'Lower'};{'Upper'};palette'];
legend([h0,h1,h2,k],legend_array,'Location','eastoutside');
ylim([0 1]);
 

