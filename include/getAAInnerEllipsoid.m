% FUNCTION TO GET THE MAXIMUM VOLUME AXIS-ALIGNED ELLIPSOID
% GIVEN THE POLYTOPE

function [E,x0] = getAAInnerEllipsoid(P)

% GET VERTICES OF POLYTOPE
V = P.V;

% FIND LIMITS IN EACH DIMENSION
nvert = size(V,1);
ndim = size(V,2);
min_vec = zeros(1,ndim);
max_vec = zeros(1,ndim);
for i = 1:ndim
  min_vec(i) = min(V(:,i));
  max_vec(i) = max(V(:,i));
end

% SCALE VERTICES
S_v = max_vec - min_vec;
for i = 1:nvert
  V(i,:) = V(i,:)./(S_v + eps);
end
P_scaled = Polyhedron(V);

% GET MAXIMUM HYPERSPHERE
data = P_scaled.chebyCenter();
x_ball = data.x; x_ball = x_ball';
R_ball = data.r;

% RESCALE AND RETURN 
x0 = S_v.*x_ball;
E = diag(R_ball*S_v);

end