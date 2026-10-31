%task5
global trace_sse trace_lik;
trace_sse = [];
trace_lik = [];

T = (0:0.001:0.010)';
R2s_true = 70;
r = [0.087,0.694,0.128,0.004,0.039,0.048]';
fF_m = [-0.498,-0.447,-0.345,-0.261,-0.063,0.064]*1e3;  % Hz
test_sets = [1.0, 0.1;
             0.1, 1.0;
             1.0, 0.3;
             0.3, 1.0];

figure('Position', [100 100 1200 1000]);
for set_idx = 1:4
    trace_sse = [];
    trace_lik = [];
    best_sse_cost = Inf;
    best_lik_cost = Inf;
    a_w_true = test_sets(set_idx,1);
    a_f_true = test_sets(set_idx,2);
    fatSig = exp(1i*2*pi*(T*fF_m)) * r;
    S_true = (a_w_true + a_f_true*fatSig) .* exp(-T*R2s_true);
    sigma_est = max(abs(S_true))/70;
    noisM = abs(S_true + sigma_est*(randn(size(S_true)) + 1i*randn(size(S_true))));
    sse_obj = @(p) sum((noisM - abs((p(1) + p(2)*fatSig).*exp(-T*p(3)))).^2);
    lik_obj = @(p) rician_nll(p, noisM, fatSig, T, sigma_est);
    nw = 100; nf = 100;
    wgrid = linspace(0,2,nw);
    fgrid = linspace(0,2,nf);
    [W,F] = meshgrid(wgrid,fgrid);
    SSEsurf = zeros(nf,nw);
    LIKsurf = zeros(nf,nw);

    fprintf('Computing surfaces for test set %d\n', set_idx);
    for i = 1:nf
        for j = 1:nw
            try
                SSEsurf(i,j) = sse_obj([W(i,j), F(i,j), R2s_true]);
                LIKsurf(i,j) = lik_obj([W(i,j), F(i,j), R2s_true]);
            catch
                SSEsurf(i,j) = Inf;
                LIKsurf(i,j) = Inf;
            end
        end
        if mod(i,10) == 0
            fprintf('Row %d of %d\n', i, nf);
        end
    end

    SSEsurf(~isfinite(SSEsurf)) = max(SSEsurf(isfinite(SSEsurf)));
    LIKsurf(~isfinite(LIKsurf)) = max(LIKsurf(isfinite(LIKsurf)));

    opts_sse = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 1000,'MaxIterations', 500,'OptimalityTolerance', 1e-6,'StepTolerance', 1e-6,'OutputFcn', @outputTrackerSSE);
    opts_lik = optimoptions('fmincon','Display', 'off','MaxFunctionEvaluations', 1000,'MaxIterations', 500,'OptimalityTolerance', 1e-6,'StepTolerance', 1e-6,'OutputFcn', @outputTrackerLIK);


    lb = [0, 0, 50];
    ub = [1, 1, 90];

    if set_idx == 1 || set_idx == 3  % Water
        p_sse = inf;
        p_lik = inf;
        ini_p = [
            a_w_true*0.8, a_f_true*1.2, R2s_true;
            a_w_true*0.9, a_f_true*1.1, R2s_true*0.95;
            a_w_true*0.7, a_f_true*1.3, R2s_true*1.05
        ];
        for i = 1:size(ini_p, 1)
            [p_sse_i, sse_i] = fmincon(sse_obj, ini_p(i,:), [], [], [], [], lb, ub, [], opts_sse);
            if sse_i < best_sse_cost
                best_sse_cost = sse_i;
                p_sse = p_sse_i;
            end
            [p_lik_i, lik_i] = fmincon(lik_obj, ini_p(i,:), [], [], [], [], lb, ub, [], opts_lik);
            if lik_i < best_lik_cost
                best_lik_cost = lik_i;
                p_lik = p_lik_i;
            end
        end

    else  % Fat
        p_sse = inf;
        p_lik = inf;
        ini_p = [
            a_w_true*1.2, a_f_true*0.8, R2s_true;
            a_w_true*1.1, a_f_true*0.9, R2s_true*0.95;
            a_w_true*1.3, a_f_true*0.7, R2s_true*1.05
        ];
        for i = 1:size(ini_p, 1)
            [p_sse_i, sse_i] = fmincon(sse_obj, ini_p(i,:), [], [], [], [], lb, ub, [], opts_sse);
            if sse_i < best_sse_cost
                best_sse_cost = sse_i;
                p_sse = p_sse_i;
            end
            [p_lik_i, lik_i] = fmincon(lik_obj, ini_p(i,:), [], [], [], [], lb, ub, [], opts_lik);
            if lik_i < best_lik_cost
                best_lik_cost = lik_i;
                p_lik = p_lik_i;
            end
        end
    end

    subplot(2,4,set_idx);
    contourf(W, F, log10(SSEsurf), 50, 'LineColor', 'none');
    colormap(parula);
    hold on;
    plot(trace_sse(:,1), trace_sse(:,2), 'r.-', 'LineWidth', 1.5);
    plot(a_w_true, a_f_true, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(p_sse(1), p_sse(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    xlabel('\a_w'); ylabel('\a_f');
    title(sprintf('Test Set %d: log_{10}(SSE) Surface\nρw=%.1f, ρf=%.1f', set_idx, a_w_true, a_f_true));
    colorbar;
    axis square;
    legend('SSE Surface', 'Optimization Path', 'True Parameters', 'SSE Solution', 'Location', 'best');

    subplot(2,4,set_idx+4);
    contourf(W, F, LIKsurf, 50, 'LineColor', 'none');
    colormap(parula);
    hold on;
    plot(trace_lik(:,1), trace_lik(:,2), 'r.-', 'LineWidth', 1.5);
    plot(a_w_true, a_f_true, 'k*', 'MarkerSize', 10, 'LineWidth', 2);
    plot(p_lik(1), p_lik(2), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    xlabel('\a_w'); ylabel('\a_f');
    title(sprintf('Test Set %d: Negative Log-Likelihood Surface\nρw=%.1f, ρf=%.1f', set_idx, a_w_true, a_f_true));
    colorbar;
    axis square;
    legend('Likelihood Surface', 'Optimization Path', 'True Parameters', 'Likelihood Solution', 'Location', 'best');
end

sgtitle('Task 5: SSE and Likelihood Surfaces with Optimization Paths for Test Sets');
saveas(gcf, 'fig/task5_surfaces.png');


function stop = outputTrackerSSE(x, optimValues, state)
    global trace_sse
    if ~strcmp(state, 'done')
        trace_sse(end+1,:) = x(1:2)';
    end
    stop = false;
end

function stop = outputTrackerLIK(x, optimValues, state)
    global trace_lik
    if ~strcmp(state, 'done')
        trace_lik(end+1,:) = x(1:2)';
    end
    stop = false;
end

function nll = rician_nll(p, noisM, fatSig, T, sigma)
    signal_model = abs((p(1) + p(2)*fatSig).*exp(-T*p(3)));
    signal_model = max(signal_model, 1e-10);
    z = noisM.*signal_model/(sigma^2);

    nll_terms = log(sigma^2) + (noisM.^2 + signal_model.^2)/(2*sigma^2);
    % Scaled besseli removes exp(abs(z)); restore it in log space.
    log_I0 = log(besseli(0, z, 1)) + abs(z);
    nll_terms = nll_terms - log_I0;
    nll = sum(nll_terms);
    if any(~isfinite(nll_terms)) || ~isfinite(nll)
        nll = 1e10;
    end
end

%% Task 6
sigma_est = max(std(diff(noisM))/sqrt(2), eps);
fprintf('Estimated noise level: %.6f\n', sigma_est);
rician_objective = @(p) rician_nll(p, noisM, fatSig, T, sigma_est);
op_r = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',2000,'MaxIterations',1000,'OptimalityTolerance',1e-6,'StepTolerance',1e-6,'FunctionTolerance',1e-6);
Nstarts = 50;
all_solutions = zeros(Nstarts,3);
all_nlls = inf(Nstarts,1);
bestNLL = Inf;
bestP_ric = [];


w_starts = [0.2, 0.5, 0.8, 1];
f_starts = [0.1, 0.3, 0.5, 1];
r2_starts = [50, 70, 90];

start_idx = 1;
for w = w_starts
    for f = f_starts
        for r2 = r2_starts
            if start_idx <= Nstarts
                p0 = [w, f, r2] .* (1 + 0.1*randn(1,3));
                p0 = max(min(p0, ub), lb);

                try
                    [p, nll] = fmincon(rician_objective, p0, [],[],[],[], lb, ub, [], op_r);
                    all_solutions(start_idx,:) = p;
                    all_nlls(start_idx) = nll;
                    if nll < bestNLL && all(isfinite(p))
                        bestNLL = nll;
                        bestP_ric = p;
                    end
                catch
                    continue;
                end
                start_idx = start_idx + 1;
            end
        end
    end
end


if isempty(bestP_ric)
    error('No finite Rician fit was obtained.');
end
success_rate = 100 * sum(abs(all_nlls - bestNLL) < 0.05*max(abs(bestNLL), eps)) / sum(isfinite(all_nlls));

% Plot results
figure('Position', [100 100 1000 400]);
subplot(1,2,1)
plot(T*1e3, noisM, 'ko', 'MarkerSize', 6); hold on;
plot(T*1e3, abs(S_true), 'k-', 'LineWidth', 2);
plot(T*1e3, abs((p_sse(1) + p_sse(2)*fatSig).*exp(-T*p_sse(3))), 'b--', 'LineWidth', 1.5);
plot(T*1e3, abs((bestP_ric(1) + bestP_ric(2)*fatSig).*exp(-T*bestP_ric(3))), 'r:', 'LineWidth', 1.5);
xlabel('Time (ms)'); ylabel('|Signal|');
legend('Data', 'True', 'SSE Fit', 'Rician Fit', 'Location', 'best');
title(sprintf('Signal Fits (σ=%.3f)', sigma_est));
grid on;

subplot(1,2,2)
valid_ix = isfinite(all_nlls);
scatter(all_solutions(valid_ix,1), all_solutions(valid_ix,2), 50,all_nlls(valid_ix)-min(all_nlls(valid_ix)), 'filled'); hold on;
plot(a_w_true, a_f_true, 'k*', 'MarkerSize', 15, 'LineWidth', 2);
plot(p_sse(1), p_sse(2), 'bs', 'MarkerSize', 10, 'LineWidth', 2);
plot(bestP_ric(1), bestP_ric(2), 'rd', 'MarkerSize', 10, 'LineWidth', 2);
colorbar; colormap(parula);
xlabel('\a_w'); ylabel('\a_f');
legend('All starts', 'True', 'SSE', 'Rician', 'Location', 'best');
title(sprintf('Parameter Space (Success: %.1f%%)', success_rate));
grid on;

sgtitle('Task 6: Rician Fitting Results');
saveas(gcf, 'fig/task6_rician.png');
