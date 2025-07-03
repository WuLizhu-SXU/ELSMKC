function WP = updateWPOPLFMVC(HP,Y,C)

k = size(HP{1}, 2);   
numker = length(HP); 
WP = cell(1, numker);

Hstar = Y * C;
for p = 1 : numker
    Tp = HP{p}' * Hstar; 
    [Up, Sp, Vp] = svd(Tp, 'econ');
    WP{p} = Up * Vp';     
end

end
