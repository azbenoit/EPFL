function execute_trace_estimation()
    % Main function to run trace estimation using GH and BGH methods
    decay_params = [0.5, 1, 1.5, 2];
    num_points = 100;
    max_exp = 3;
    gh_errors = zeros(length(decay_params), num_points);
    bgh_errors = zeros(length(decay_params), num_points);
    
    % Loop over each decay parameter to evaluate performance
    for param_idx = 1:length(decay_params)
        decay_param = decay_params(param_idx);
        gaussian_mat = randn(1000, 1000);
        [orth_matrix, ~] = qr(gaussian_mat);
        matrix_L = orth_matrix' * diag((1:1000).^-decay_param) * orth_matrix;
        clear orth_matrix;
        exact_trace = trace(matrix_L);
        
        % Loop over each number of random vectors
        for point_idx = 1:num_points
            num_random_vectors = round(logspace(1, max_exp, num_points));
            trace_gh = estimate_trace_gh(matrix_L, num_random_vectors(point_idx));
            trace_bgh = estimate_trace_bgh(matrix_L, num_random_vectors(point_idx));
            gh_errors(param_idx, point_idx) = abs(exact_trace - trace_gh) / abs(exact_trace);
            bgh_errors(param_idx, point_idx) = abs(exact_trace - trace_bgh) / abs(exact_trace);
        end
        clear matrix_L;
        
        % Plot individual results for each decay parameter
        figure;
        loglog(logspace(1, max_exp, num_points), gh_errors(param_idx, :), '-o', 'DisplayName', sprintf('GH, decay=%g', decay_param));
        hold on;
        loglog(logspace(1, max_exp, num_points), bgh_errors(param_idx, :), '-o', 'DisplayName', sprintf('BGH, decay=%g', decay_param));
        legend('Location', 'southwest');
        xlabel('Number of Samples');
        ylabel('Relative Error');
        title(sprintf('Trace Estimation Comparison for Decay %g', decay_param));
        hold off;
        
        % Save individual plot
        saveas(gcf, sprintf('trace_estimation_decay_%g.png', decay_param));
    end
    
    % Combined plot for all decay parameters
    figure;
    for param_idx = 1:length(decay_params)
        decay_param = decay_params(param_idx);
        loglog(logspace(1, max_exp, num_points), gh_errors(param_idx, :), 'DisplayName', sprintf('GH, decay=%g', decay_param));
        hold on;
        loglog(logspace(1, max_exp, num_points), bgh_errors(param_idx, :), 'DisplayName', sprintf('BGH, decay=%g', decay_param));
    end
    legend('Location', 'southwest');
    xlabel('Number of Samples');
    ylabel('Relative Error');
    title('Combined Trace Estimation for All Decay Values');
    hold off;
    
    % Save combined plot
    saveas(gcf, 'trace_estimation_combined.png');
end

function trace_est = estimate_trace_gh(mat, num_samples)
    % Girard-Hutchinson trace estimator
    trace_accum = 0;
    for sample_idx = 1:num_samples
        rand_vec = randn(size(mat, 1), 1);
        trace_accum = trace_accum + rand_vec' * mat * rand_vec;
    end
    trace_est = trace_accum / num_samples;
end

function trace_est = estimate_trace_bgh(mat, num_samples)
    % Improved Girard-Hutchinson trace estimator
    rank = num_samples;
    proj_rank = rank + 1;
    add_rank = 2 * (rank + 1);
    omega = randn(size(mat, 1), rank + proj_rank);
    psi = randn(size(mat, 1), rank + proj_rank + add_rank);
    [ortho_basis, ~] = qr(mat * omega, 0);
    clear omega;
    proj_mat = ortho_basis * ((psi' * ortho_basis) \ psi');
    clear psi ortho_basis;
    approx_mat = proj_mat * mat;
    remainder = mat - approx_mat;
    trace_approx = trace(approx_mat);
    clear approx_mat;
    
    trace_accum = 0;
    for sample_idx = 1:num_samples
        rand_vec = randn(size(mat, 1), 1);
        trace_accum = trace_accum + rand_vec' * remainder * rand_vec;
    end
    trace_est = (trace_accum / num_samples) + trace_approx;
end
