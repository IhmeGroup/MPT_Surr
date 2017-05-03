function [x0_out,diff_out] = getOuterCuboid(A,b)

m = size(A,1);
n = size(A,2);
l_out = zeros(1,n); u_out = zeros(1,n);

for i = 1:n
   cvx_begin quiet
     variables x_out(n);
     minimize x_out(i)
     subject to
       A*x_out <= b
   cvx_end
   
   l_out(i) = x_out(i);
   
   cvx_begin quiet
     variables x_out(n);
     maximize x_out(i)
     subject to
       A*x_out <= b
   cvx_end
   
   u_out(i) = x_out(i);
    
end

x0_out = 0.5*(l_out+u_out);
diff_out = 0.5*(u_out-l_out);

end