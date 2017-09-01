function [palette,palette_label,palette_tsi,palette_dcn,exp_comp,target_mw,target_hc] = ...
    set_target_details(surr)


% PALETTE COMPOUNDS
switch surr
    case 'dooley'
        palette = {'NPBENZ','NC12H26','IC8H18','TMBENZ'};
    otherwise
        palette = {'MCYC6','C7H8','C6H6','IC8H18','NC12H26'};
end

% LABELS
switch surr
    case 'dooley'
        palette_label = {'n-Propylbenzene', 'n-Dodecane', 'iso-Octane', 'Trimethylbenzene'};
    otherwise
        palette_label = {'Methyl-Cyclohexane', 'Toluene', 'Benzene', 'iso-Octane', 'n-Dodecane'};
end

% PALETTE TSI
tsi = load_tsi_database();
palette_tsi = zeros(length(palette),1);
for i = 1:length(palette)
    palette_tsi(i) =tsi.(palette{i});
end

% EXPERIMENTAL COMPOSITION (MOLE FRACTIONS)
switch surr
    case 'dooley'
        exp_comp = [0.228, 0.404, 0.295, 0.073];
    case 'violi'
        exp_comp = [0.1, 0.1, 0.01, 0.055, 0.735];
    case 'hanson_a'
        exp_comp = [0.1, 0.1, 0.01, 0.25, 0.54];
    case 'hanson_b'
        exp_comp = [0.1, 0.295, 0.01, 0.055, 0.54];
end

% PURE COMPONENT IDT
switch surr
    % TAKEN FROM M.J. Murphy, J.D. Taylor, R.L. McCormick, Compendium of Experimental Cetane
    % Number Data, National Renewable Energy Laboratory, 2004.
    case 'dooley'
        palette_dcn = [28.2,78,17,21.8]';
    otherwise
        palette_dcn = [20,9,23,14,78]';
end

% % TARGET PROPERTIES
switch surr
    case 'dooley'
        target_mw = 138.7;
    case 'test'
        target_mw = 138.7;
    case 'violi'
        target_mw = 151.3;
    case 'hanson_a'
        target_mw = 140.3;
    case 'hanson_b'
        target_mw = 136.0;
end

switch surr
    case 'dooley'
        target_hc = 1.96;
    case 'test'
        target_hc = 1.96;
    case 'violi'
        target_hc = 2.08;
    case 'hanson_a'
        target_hc = 2.09;
    case 'hanson_b'
        target_hc = 1.93;
end

end