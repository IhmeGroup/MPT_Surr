function draw_both_2d(P,Ell,Ell_AA,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label,save_output,surr)

palette_size = length(palette_label);
tuples = gen_tuples(1:palette_size-1,2);

for i = 1:size(tuples,1)

figure('Position', [100, 100, 1049, 895]);    
    
% PROJECTION BASIS
idx_1 = tuples(i,1);
idx_2 = tuples(i,2);

% ORTHOGONAL MATRIX
basis = zeros(palette_size-1,2);
basis(idx_1,1) = 1;
basis(idx_2,2) = 1;  


dims = [idx_1 idx_2];
P_proj = P.projection(dims);

% ELLIPSOIDS
Ell_AA_Proj = projection(Ell_AA,basis);
Ell_Proj = projection(Ell,basis);
opts=[];
opts.fill = [0 0];
opts.width = [4 2];
opts.style = [':','-'];
opts.color = [0 0 0;0 0 0];
plot(Ell_AA_Proj,Ell_Proj,opts);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',40,'FontName','Times New Roman');

hold on;

a = P_proj.plot();
set(a,'LineWidth',2)
set(a,'EdgeColor',[1 0 0],'FaceColor',[1 1 1]);
alpha(a,0);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',40,'FontName','Times New Roman');
grid off;

% DRAW HYPERCUBOIDS
inner = rectangle('Position',[x0_R(idx_1)-diff_R(idx_1) x0_R(idx_2)-diff_R(idx_2) 2*diff_R(idx_1) 2*diff_R(idx_2)],'LineWidth',2);
outer = rectangle('Position',[x0_out(idx_1)-diff_out(idx_1) x0_out(idx_2)-diff_out(idx_2) 2*diff_out(idx_1) 2*diff_out(idx_2)],'LineWidth',2,'EdgeColor',[0 0 1]);

h = plot(exp_comp(idx_1),exp_comp(idx_2),'^','MarkerSize',40);
set(h,'MarkerEdgeColor',[1 0 0 ],'MarkerFaceColor',[1 0 0 ]);

% SET LABELS
xlabel(palette_label{idx_1});
ylabel(palette_label{idx_2});

if (save_output)
file_name = strcat('./Plots/WonThreshold/',surr,'/LC',num2str(i));
savefig(file_name);
disp(strcat('Written file: ',file_name));

file_name = strcat('./Images/WonThreshold/',surr,'/LC',num2str(i),'.pdf');
export_fig(file_name);
disp(strcat('Written file: ',file_name));
end

end

end