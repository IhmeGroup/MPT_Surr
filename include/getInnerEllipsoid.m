function [Ein,cin,yin] = getInnerEllipsoid_cvx(P)

n = P.Dim;
m = size(P.A,1);

cvx_begin
    variable B(n,n) symmetric
    variable d(n)
    dual variable y{m}
    maximize( det_rootn( B ) )
    subject to
       for i = 1:m
           y{i} : norm( B*P.A(i,:)', 2 ) + P.A(i,:)*d <= P.b(i);
       end
cvx_end

Ein = double(B);
cin = double(d);
yin = double(cell2mat(y));

end