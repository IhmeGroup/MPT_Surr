function draw_ellipsoids_2d(P,Ell,Ell_AA,exp_comp,palette_label,save_output,surr)

palette_size = length(palette_label);
tuples = gen_tuples(1:palette_size-1,2);

for i = 1:size(tuples,1)

% PROJECTION BASIS
idx_1 = tuples(i,1);
idx_2 = tuples(i,2);

% ORTHOGONAL MATRIX
basis = zeros(palette_size-1,2);
basis(idx_1,1) = 1;
basis(idx_2,2) = 1;  

% CREATE PROJECTIONS
Ell_AA_Proj = projection(Ell_AA,basis);
Ell_Proj = projection(Ell,basis);
    
figure('Position', [100, 100, 1049, 895]);    
hold on;

dims = [idx_1 idx_2];
P_proj = P.projection(dims);

a = P_proj.plot();
set(a,'LineWidth',4)
set(a,'EdgeColor',[0 0 0],'FaceColor',[1 1 1]);
alpha(a,0);
grid off;

opts=[];
opts.fill = [0 0];
opts.width = [1 4];
opts.style = [':','-'];
opts.color = [0 0 0;0 0 0];
plot(Ell_AA_Proj,Ell_Proj,opts);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',40,'FontName','Times New Roman');

h = plot(exp_comp(idx_1),exp_comp(idx_2),'^','MarkerSize',40);
set(h,'MarkerEdgeColor',[1 0 0 ],'MarkerFaceColor',[1 0 0 ]);

% SET LABELS
xlabel(palette_label{idx_1});
ylabel(palette_label{idx_2});

if (save_output)
file_name = strcat('./Plots/',surr,'_hyperellipsoid/LC',num2str(i));
savefig(file_name);
disp(strcat('Written file: ',file_name));

file_name = strcat('./Images/',surr,'_hyperellipsoid/LC',num2str(i),'.pdf');
export_fig(file_name);
disp(strcat('Written file: ',file_name));
end

end

end