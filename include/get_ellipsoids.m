function [Ell_AA,E_AA,Ell,E,x0_AA,x0,y0] = get_ellipsoids(P)

% GET AXIS-ALIGNED INNER ELLIPSOID
[E_AA,x0_AA] = getAAInnerEllipsoid(P);

% GET INNER ELLIPSOID
[E,x0,y0] = getInnerEllipsoid(P);x0 = x0';y0=[y0(1) y0(3)];y0 = y0/norm(y0);

% PROJECTIONS
% ELLIPSOIDS DESCRIBED BY (x-q)Q^-1(x-q)
% Also E'E = Q gives Ex+F form
Ell_AA = ellipsoid(x0_AA',E_AA*E_AA'); Ell = ellipsoid(x0',E*E');
disp('Dual variables');
disp(y0);

end