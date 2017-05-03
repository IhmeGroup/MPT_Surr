function [tuples] = gen_tuples(vec,k)

% ASSUMING vec ENTRIES ARE ALWAYS DISTINCT
palette_size = length(vec);
tuples = [];
vec = sort(vec);

% GENERATE PAIRWISE INCREASING FOR 2
if (k == 1)
  for i = 1:palette_size
      tuples = [tuples;vec(i)];
  end

% PICK AN ELEMENT AND RECURSIVELY APPLY
else
  for p = 1:palette_size
    buf = gen_tuples([vec(p+1:end)],k-1);
    if (length(buf) ~= 0)
      tuples = [tuples; [vec(p)*ones(size(buf,1),1) buf]];
    else
      tuples = [tuples; vec(p)*ones(size(buf,1),1)];  
  end
end

end