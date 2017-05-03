function [vol] = sp_vol(n)

if (mod(n,2)==0)
   vol = (pi^(n/2))/factorial(n/2);
else
   n_h = (n-1)/2;
   vol = (2^n)*(pi^n_h)*factorial(n_h)/factorial(n);
end

end