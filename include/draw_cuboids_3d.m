function draw_cuboids_3d(P,x0_R,diff_R,x0_out,diff_out,exp_comp,palette_label)

palette_size = length(palette_label);
tuples = gen_tuples(1:palette_size-1,3);

for i = 1:size(tuples,1)

figure('Position', [100, 100, 1049, 895]);    
    
idx_1 = tuples(i,1);
idx_2 = tuples(i,2);
idx_3 = tuples(i,3);

dims = [idx_1 idx_2 idx_3];
P_proj = P.projection(dims);
    
% POLYTOPE

a = plot(P_proj);
set(a,'LineWidth',4)
alpha(a,0);
set(gcf,'Color',[1 1 1]);
set(gca,'FontSize',32,'FontName','Times New Roman');
grid off;

hold on;

% CREATE TRUNCATED COORDINATES
diff_RT = [diff_R(idx_1) diff_R(idx_2) diff_R(idx_3)];
diff_outT = [diff_out(idx_1) diff_out(idx_2) diff_out(idx_3)];
x0_RT = [x0_R(idx_1) x0_R(idx_2) x0_R(idx_3)];
x0_outT = [x0_out(idx_1) x0_out(idx_2) x0_out(idx_3)];

% INNER CUBOID
[h_in,~,~] = DrawCuboid(2*diff_RT',x0_RT',[0;0;0],'y',0.5);
set(h_in,'LineWidth',2)

% OUTER CUBOID
[h_out,~,~] = DrawCuboid(2*diff_outT',x0_outT',[0;0;0],'b',0);
set(h_out,'EdgeColor',[1 0 0],'LineWidth',2);

% SET LABELS
xlabel(palette_label{idx_1});
ylabel(palette_label{idx_2});
zlabel(palette_label{idx_3});

% EXPERIMENTAL
h = plot3(exp_comp(idx_1),exp_comp(idx_2),exp_comp(idx_3),'s','MarkerSize',40);
set(h,'MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[0 0 1]);

end

end