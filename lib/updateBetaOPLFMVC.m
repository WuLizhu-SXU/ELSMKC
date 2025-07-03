function beta = updateBetaOPLFMVC(HP, WP, Y, C)
numker = length(WP);
HHPWP = zeros(numker, 1);
Hstar = Y * C;
for p = 1:numker
    if isempty(HP{p}) || isempty(WP{p})
        error('核矩阵 HP{%d} 或权重矩阵 WP{%d} 未定义或为空', p, p);
    end
    if size(WP{p}, 1) ~= size(HP{p}, 2)
        error('WP{%d} 和 HP{%d} 的维度不匹配', p, p);
    end
    HHPWP(p) = trace(Hstar' * (HP{p} * WP{p}));
end
normHHPWP = norm(HHPWP);
if normHHPWP == 0
    error('HHPWP 的范数为零，无法归一化');
end
beta = HHPWP ./ normHHPWP;
beta(beta < eps) = 0;
normBeta = norm(beta);
if normBeta == 0
    error('剔除小值后 beta 的范数为零，无法归一化');
end
beta = beta ./ normBeta;
end
