%% SOLVE MPT PROGRAM FOR MW AND HC RATIO
% AUTHOR : G. PAVAN BHARADWAJ

clear all;
clc; close all; yalmip('clear');cvx_clear;

%% INPUTS

% SURROGATE NAME
surr = 'hanson_a';

% TARGET PROPERTIES FLAG
useCompositionOnly = true;

% SAVE OUTPUT
save_output = 0;

%% LOAD MECHANISM
adddir('/Users/gpavanb/Desktop/Academics/Stanford/Ihme_Research/Surrogates/Mechanisms/POLIMI_TOT');
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
  target_tsi = exp_comp*palette_tsi;
end

%% MPT PROBLEM PARAMETERS

% ERROR THRESHOLDS
MW_err = 25.0;
TSI_err = 4.0;

% MPT VARIABLES 
x = sdpvar(palette_size, 1);
err = sdpvar(1, 2);

% CONSTRUCT CONSTRAINT MATRIX
b = [target_mw + err(1); -target_mw + err(1);target_tsi + err(2); -target_tsi + err(2)];
A = [MW_Vec';-MW_Vec';palette_tsi';-palette_tsi'];

% CONSTRAINTS
C = [A*x <= b, x >= 0, sum(x)==1, 0 <= err(1) <= MW_err, 0 <= err(2) <= TSI_err];

%% SOLVE MPT PROBLEM
x_ax = linspace(0,MW_err,20);
y_ax = linspace(0,TSI_err,20);

colors = distinguishable_colors(palette_size,'g');
figure('Position', [100, 100, 1000, 895]);

for i = 1:palette_size

    % OBJECTIVE FUNCTION
    sel = zeros(1,palette_size);
    sel(i) = 1;
    J = sel*x;
    
    % SOLVE MINIMUM PROBLEM
    plp = Opt(C, J, err, x);
    solution = plp.solve();
    
    for j = 1:length(x_ax)
        for k = 1:length(y_ax)
            Zl(j,k) = solution.xopt.feval([x_ax(j);y_ax(k)],'obj');
        end
    end
    
    sel = zeros(1,palette_size);
    sel(i) = -1;
    J = sel*x;
    
    % SOLVE MAXIMUM PROBLEM
    plp = Opt(C, J, err, x);
    solution = plp.solve();
    
    for j = 1:length(x_ax)
        for k = 1:length(y_ax)
            Zu(j,k) = -solution.xopt.feval([x_ax(j);y_ax(k)],'obj');
        end
    end
    
    % 3D PLOT
    
    hold all;
    [X,Y] = meshgrid(x_ax,y_ax);
    hsl = surf(X,Y,Zl);
    set(hsl,'FaceColor',colors(i,:),'EdgeColor','none');
    hsu = surf(X,Y,Zu);
    set(hsu,'FaceColor',colors(i,:),'EdgeColor','none');

end


set(gcf,'Color',[1 1 1]);
set(gca,'FontName','Times New Roman','FontSize',32);
xlabel('\epsilon_{MW}','Interpreter','Tex');
ylabel('\epsilon_{TSI}','Interpreter','Tex');
zlabel('Mole Fraction');

% EXPERIMENTAL COMPOSITION
for j = 1:palette_size
  h = plot3(0,0,exp_comp(j),'^','MarkerSize',20);
  set(h,'MarkerEdgeColor',colors(j,:),'MarkerFaceColor',colors(j,:));
end

h0 = plot3(-1,-1,0,'^','MarkerSize',20,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
legend_array = [{'Experimental'};palette'];

% UPPER BOUNDS
for j = 1:palette_size
    legs(j) = plot3(-1,-1,0,'-','LineWidth',2,'Color',colors(j,:));
end

legend([h0,legs],legend_array,'Location','eastoutside');
xlim([0 MW_err]);
ylim([0 TSI_err]);
zlim([0 1]);

if (save_output)
    file_name = strcat('./Plots/WonThreshold/mpt/',surr);
    savefig(file_name);
    disp(strcat('Written file: ',file_name));
    
    file_name = strcat('./Images/WonThreshold/mpt/',surr,'.pdf');
    export_fig(file_name);
    disp(strcat('Written file: ',file_name));
end

