function disp_ellipse_stats(P,x0_AA,E_AA,x0,E)

x0_AA = [x0_AA,1-sum(x0_AA)];
x0 = [x0,1-sum(x0)];

palette_size = length(x0_AA);

% DISPLAY FINAL PROPERTIES

disp('AA center');
disp(num2str(x0_AA));
if (sum(x0_AA(1:end-1))>1)
  disp('Do not use center. Other feasible points exist');
end
disp('Intervals');
disp(num2str(diag(E_AA)));
disp(' ');
disp('MV center');
disp(num2str(x0));
if (sum(x0(1:end-1))>1)
  disp('Do not use center. Other feasible points exist');
end
disp('Matrix');
disp(num2str(E));
disp(' ');

% PRINT VOLUME STATISTICS
poly_vol = P.volume();
AA_vol = sp_vol(palette_size-1)*det(E_AA);
MV_vol = sp_vol(palette_size-1)*det(E);

disp('Axis-Aligned Ellipse');
disp(strcat('Packing fraction: ',num2str(AA_vol/poly_vol)));
disp('Maximum Volume Ellipse');
disp(strcat('Packing fraction: ',num2str(MV_vol/poly_vol)));

end