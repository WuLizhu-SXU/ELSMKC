clear;
clc;

data_path = fullfile(pwd, 'data');
addpath(data_path);
lib_path = fullfile(pwd, 'lib');
addpath(lib_path);
dirop = dir(fullfile(data_path, '*.mat'));
datasetCandi = {dirop.name};
exp_n = 'ELSMKC';


for i1 = 1 : length(datasetCandi)
    data_name = datasetCandi{i1}(1:end-4);
    dir_name = fullfile(pwd, exp_n, data_name);

    try
        if ~exist(dir_name, 'dir')
            mkdir(dir_name);
        end
        prefix_mdcs = dir_name;
    catch
        disp(['Create dir: ', dir_name, ' failed, check authorization']);
    end

    clear X y Y;

    load(data_name);
    if exist('y', 'var')
        Y = y;
    end
    if size(X, 1) ~= size(Y, 1)
        Y = Y';
    end
    assert(size(X, 1) == size(Y, 1));
    nSmp = size(X, 1);
    nCluster = length(unique(Y));

    fname2 = fullfile(prefix_mdcs, [data_name, '_ELSMKC.mat']);
    if ~exist(fname2, 'file')
        Xs = {X};
        Ks = Xs_to_Ks_12k(Xs);
        Ks = Ks{1, 1};
        nKernel = size(Ks, 3);
        Ks_cell = cell(1, nKernel);
        for iKernel = 1:nKernel
            Ks_cell{iKernel} = Ks(:, :, iKernel);
        end
        clear Ks;
       lambdas = 2.^(-5:1:5);
        nRepeat = 10;
        nMeasures = 13;
        seed = 2024;
        rng(seed);
        random_seeds = randi([0, 1e6], 1, nRepeat);
        original_rng_state = rng;

        nParam = length(lambdas);
        ELSMKC_result = zeros(nParam, 1, nRepeat, nMeasures);
        ELSMKC_time = zeros(nParam, 1);
     
        % 核、特征初始化只做一次
        t1_s = tic;
        LKs_cell = Ks_cell;
        [Hs, LKHs] = Ks2Hs(LKs_cell, nCluster);
        [avgH, avgH_normalized] = LKs2AvgH(LKs_cell, nCluster);
        t1 = toc(t1_s);

        for iParam = 1:length(lambdas)
            lambda = lambdas(iParam);
            t2 = 0;

            for iRepeat = 1:nRepeat
                disp(['Lambda ', num2str(lambda), ' Repeat ', num2str(iRepeat)]);
                rng(original_rng_state);
                rng(random_seeds(iRepeat));
                t2_s = tic;
                label0= litekmeans(avgH_normalized, nCluster, 'MaxIter', 100, 'Replicates', 10);
                Y0 = ind2vec(label0')';
                [label, Ws, C, alpha, beta, objHistory] = ELSMKC(LKs_cell, Hs, LKHs,lambda, Y0);
                t2 = t2 + toc(t2_s);
                result_i = my_eval_y(label, Y)';
               ELSMKC_result(iParam, 1, iRepeat, :) = result_i;
            end

            tt = t1 + t2 / nRepeat;
            ts = [t1, t2];
            ELSMKC_time(iParam) = tt;

            fname3 = fullfile(prefix_mdcs, [data_name, '_lambda_', num2str(iParam), '.mat']);
            save(fname3, 'ELSMKC_result', 'ts');
        end
      
        a1 = sum(ELSMKC_result, 3);
        a2 = reshape(a1, nParam, nMeasures);
        ELSMKC_result_grid = [a2 / nRepeat, ELSMKC_time];
        ELSMKC_result_summary = [max(ELSMKC_result_grid, [], 1), mean(ELSMKC_time)];

        save(fname2, 'ELSMKC_result', 'ELSMKC_result_grid', 'ELSMKC_time', 'ELSMKC_result_summary');
        disp([data_name(1:end-4), ' has been completed!']);
    end
end
rmpath(data_path);
rmpath(lib_path);
