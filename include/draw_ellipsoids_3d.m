function draw_ellipsoids_3d(P,Ell,Ell_AA,exp_comp,palette_label)

palette_size = length(palette_label);
tuples = gen_tuples(1:palette_size-1,3);

for i = 1:size(tuples,1)

% PROJECTION BASIS
idx_1 = tuples(i,1);
idx_2 = tuples(i,2);
idx_3 = tuples(i,3);

% ORTHOGONAL MATRIX
basis = zeros(palette_size-1,3);
basis(idx_1,1) = 1;
basis(idx_2,2) = 1; 
basis(idx_3,3) = 1;

% CREATE PROJECTIONS
Ell_AA_Proj = projection(Ell_AA,basis);
Ell_Proj = projection(Ell,basis);
    
figure('Position', [100, 100, 1049, 895]);    
hold on;

dims = [idx_1 idx_2 idx_3];
P_proj = P.projection(dims);

opts=[];
opts.fill = [0 0];
opts.width = [1 4];
opts.style = ['--','--'];
opts.color = [0 0 1;1 1 0];
opts.shade = [0.4 0.4];
plot(Ell_AA_Proj,Ell_Proj,opts);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',40,'FontName','Times New Roman');

h = plot3(exp_comp(idx_1),exp_comp(idx_2),exp_comp(idx_3),'^','MarkerSize',40);
set(h,'MarkerEdgeColor',[1 0 0 ],'MarkerFaceColor',[1 0 0 ]);

a = P_proj.plot();
set(a,'LineWidth',4)
alpha(a,0);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',32,'FontName','Times New Roman');
grid off;

% SET LABELS
xlabel(palette_label{idx_1});
ylabel(palette_label{idx_2});
zlabel(palette_label{idx_3});

end
    
end