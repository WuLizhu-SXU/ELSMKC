function [LKs, Ks_avg] = Ks2LKs(Ks_cell, knn_size)

nSmp = size(Ks_cell{1}, 1);
nKernel = length(Ks_cell);

Ks_avg = Ks_cell{1};
for iKernel = 2:nKernel
    Ks_avg = Ks_avg + Ks_cell{iKernel};
end
Ks_avg = Ks_avg/nKernel;

NS = genarate_neighborhood(Ks_avg, knn_size);%% k*nSmp

Ai_sum = zeros(nSmp);
for iSmp = 1:nSmp
    Ai_sum(NS(1:knn_size, iSmp), NS(1:knn_size, iSmp)) = Ai_sum(NS(1:knn_size, iSmp), NS(1:knn_size, iSmp)) + 1;
end
Ai_sum = Ai_sum./nSmp;

LKs = cell(1, nKernel);
for iKernel=1:nKernel
    LKs{iKernel} = Ks_cell{iKernel} .* Ai_sum;
end
end

function [indx_0]  = genarate_neighborhood(KC,k)
nSmp = size(KC,1);
KC0 = KC - 10^8*eye(nSmp);
[~,indx] = sort(KC0,'descend');
indx_0 = indx(1:k,:);
end