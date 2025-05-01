function [is_fault, excluded_idx] = raim_detection(A, omc, weight, settings)
% RAIM 故障检测与排除
% 输入：
%   A       - 几何矩阵
%   omc     - 观测值残差 (obs - predicted)
%   weight  - 权重向量 (diagonal of C)
% 输出：
%   is_fault      - 是否检测到故障 (true/false)
%   excluded_idx  - 建议排除的卫星索引

    
    W = diag(weight);
    %W = diag(ones(size(A,1),1))/(settings.sigma^2);
    



    residuals = omc;
    T = residuals' * W * residuals;
    sigma = settings.sigma;
    T_threshold = 5.33 * sigma;

    if T > T_threshold
        min_T = inf;
        excluded_idx = 0;
        for i = 1:size(A, 1)
            % 排除第 i 颗卫星
            A_excluded = A; A_excluded(i, :) = [];
            W_excluded = W; W_excluded(i, :) = []; W_excluded(:, i) = [];
            omc_excluded = omc; omc_excluded(i) = [];
            
            % 检查矩阵秩并正则化
            H_excluded = A_excluded' * W_excluded * A_excluded;
            if rank(H_excluded) < 4
                continue;  % 跳过几何病态的卫星排除
            end
            x_excluded = (H_excluded + 1e-8 * eye(4)) \ (A_excluded' * W_excluded * omc_excluded);
            
            % 计算残差和检验统计量
            residuals_excluded = omc_excluded - A_excluded * x_excluded;
            T_i = residuals_excluded' * W_excluded * residuals_excluded;
            
            if T_i < min_T
                min_T = T_i;
                excluded_idx = i;
            end
        end
        is_fault = (excluded_idx > 0);  % 仅当找到有效排除时标记故障
    else
        is_fault = false;
        excluded_idx = 0;
    end
end