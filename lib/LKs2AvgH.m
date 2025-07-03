function [avgH, avgH_normalized] = LKs2AvgH(LKs, nCluster)
nKernel = length(LKs);
nSmp = size(LKs{1}, 1);
LK_avg = zeros(nSmp, nSmp);
for iKernel = 1:nKernel
    LK_avg = LK_avg + LKs{iKernel};
end
LK_avg = (LK_avg + LK_avg')/2;
opt.disp = 0;
[avgH, ~] = eigs(LK_avg, nCluster, 'LA', opt);

avgH_normalized = avgH ./ (repmat(sqrt(sum(avgH.^2, 2)), 1, nCluster) + eps);
end