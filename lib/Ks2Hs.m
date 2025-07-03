function [Hs, KHs] = Ks2Hs(Ks_cell, nCluster)
nKernel = length(Ks_cell);
Hs = cell(1, nKernel);
if nargout == 2
    KHs = cell(1, nKernel);
end
opt.disp = 0;
for iKernel = 1:nKernel
    Ki = Ks_cell{iKernel};
    Ki = (Ki+Ki')/2;
    [Hs{iKernel}, ~] = eigs(Ki, nCluster, 'LA', opt);
    if nargout == 2
        KHs{iKernel} = Ki * Hs{iKernel};
    end
end
end