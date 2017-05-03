function [Ein,cin] = getInnerEllipsoid(P)

E = sdpvar(P.Dim, P.Dim); 
c = sdpvar(P.Dim, 1); 
con = []; 
for i = 1:size(P.A, 1) 
    con = con + [ norm(E*P.A(i, :)') + P.A(i, :)*c <= P.b(i) ]; 
end 
con = con + [ E>=0 ];
opt = sdpsettings('solver','sdpt3','verbose',0);
optimize(con, -logdet(E));
Ein = double(E);
disp('');
disp(num2str(-log(det(Ein))));
disp('');
cin = double(c);

end