function [x0_R,diff_R,y] = getInnerCuboid(A,b)

A_pos = A;
A_neg = A;
m = size(A,1);
n = size(A,2);
for i=1:m
    for j=1:n
        A_pos(i,j) = max(A(i,j), 0);
        A_neg(i,j) = min(A(i,j), 0);
    end
end

cvx_begin
    variables l(n) u(n);
    dual variable y;
    maximize geo_mean(u - l);
    subject to
        y : A_pos*u + A_neg*l <= b % because the LHS is as big as it could ever be
        u >= 0
        l >= 0
cvx_end

% RETURN CENTER AND TOLERANCE RADIUS
l=l';u=u';
x0_R = 0.5*(l+u);
diff_R = 0.5*(u-l);

end