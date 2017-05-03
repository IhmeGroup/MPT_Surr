function disp_cuboid_stats(P,x0_R,diff_R,x0_out,diff_out,y)

% PRINT STATISTICS
x0_R = [x0_R,1-sum(x0_R)];
x0_out = [x0_out,1-sum(x0_out)];

% DISPLAY FINAL PROPERTIES
disp('Inner Hypercuboid center');
disp(num2str(x0_R));
if (sum(x0_R(1:end-1))>1)
  disp('Do not use center. Other feasible points exist');
end
disp('Intervals');
disp(num2str(diff_R));
disp(' ');
disp('Outer Hypercuboid center');
disp(num2str(x0_out));
if (sum(x0_out(1:end-1))>1)
  disp('Do not use center. Other feasible points exist');
end
disp('Intervals');
disp(num2str(diff_out));
disp(' ');
disp('Dual Variables');
y_fin = [y(1) y(3)];
disp(y_fin/norm(y_fin));

% PRINT VOLUME STATISTICS
poly_vol = P.volume();
rect_in_vol = prod(2*diff_R);
rect_out_vol = prod(2*diff_out);

disp(strcat('Packing fraction: ',num2str(rect_in_vol/poly_vol)));
disp(strcat('Outer to Inner: ',num2str(rect_out_vol/rect_in_vol)));

end