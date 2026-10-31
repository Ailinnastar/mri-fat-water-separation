%% Task 7

T = (0:0.001:0.010)';
SNR = 70;
fF_m = [-0.498,-0.447,-0.345,-0.261,-0.063,0.064]*1e3;
r = [0.087,0.694,0.128,0.004,0.039,0.048]';
apw_values = [0.2, 0.5, 0.8, 1.0];
apf_values = [0.1, 0.3, 0.5, 0.7];

n_w = length(apw_values);
n_f = length(apf_values);
R2s_water_errors = zeros(n_w, n_f);
R2s_fat_errors = zeros(n_w, n_f);
fit_success_rates = zeros(n_w, n_f);

R2s_water_true = 70;
R2s_fat_true = 90;

fatSignal = exp(1i*2*pi*(T*fF_m)) * r;

Nstarts_total = 20;
opts = optimoptions('fmincon','Display', 'off','MaxIterations', 1000,'OptimalityTolerance', 1e-6,'StepTolerance', 1e-6);

lb = [0.001, 0.001, 50, 50];
ub = [1, 1, 100, 100];

for i = 1:n_w
    for j = 1:n_f

        S_true = (apw_values(i) * exp(-T*R2s_water_true) + apf_values(j) * fatSignal .* exp(-T*R2s_fat_true));
        noise_level = max(abs(S_true))/SNR;
        S_noisy = S_true + noise_level*(randn(size(S_true)) + 1i*randn(size(S_true)));

        all_params = zeros(Nstarts_total, 4);
        all_costs = zeros(Nstarts_total, 1);

        for k = 1:Nstarts_total
            if rand < 0.5
                p0 = [apw_values(i)*(0.8+0.4*rand), apf_values(j)*(0.6+0.8*rand), R2s_water_true*(0.8+0.4*rand), R2s_fat_true*(0.8+0.4*rand)];
            else
                p0 = [apw_values(i)*(0.6+0.8*rand), apf_values(j)*(0.8+0.4*rand), R2s_water_true*(0.8+0.4*rand), R2s_fat_true*(0.8+0.4*rand)];
            end

            p0 = max(lb, min(ub, p0));
            objfun = @(params) sum(abs(abs(S_noisy) - abs((params(1)*exp(-T*params(3)) + params(2)*fatSignal.*exp(-T*params(4))))).^2);

            [params, cost] = fmincon(objfun, p0, [], [], [], [], lb, ub, [], opts);

            all_params(k,:) = params;
            all_costs(k,:) = cost;
        end

        [~, best_idx] = min(all_costs);
        best_p = all_params(best_idx,:);

        denom = max(abs(best_p), 1e-6);
        param_errors = abs(all_params - best_p) ./ denom;
        similar_to_best = sum(param_errors < 0.05, 2) == 4;

        R2s_water_errors(i,j) = 100*abs(best_p(3) - R2s_water_true)/R2s_water_true;
        R2s_fat_errors(i,j) = 100*abs(best_p(4) - R2s_fat_true)/R2s_fat_true;
        fit_success_rates(i,j) = 100 * sum(similar_to_best) / Nstarts_total;

    end
end

figure('Name', 'Task 7: Dual R2* Stability Analysis', 'Position', [100 100 1200 800]);

subplot(2,2,1);
imagesc(apf_values, apw_values, R2s_water_errors);
colormap(jet); colorbar;
xlabel('Fat Amplitude (ρf)');
ylabel('Water Amplitude (ρw)');
title('R2*_{water} Error (%)');
axis square;


subplot(2,2,2);
imagesc(apf_values, apw_values, R2s_fat_errors);
colormap(jet); colorbar;
xlabel('Fat Amplitude (ρf)');
ylabel('Water Amplitude (ρw)');
title('R2*_{fat} Error (%)');
axis square;

subplot(2,2,3);
imagesc(apf_values, apw_values, fit_success_rates);
colormap(jet); colorbar;
xlabel('Fat Amplitude (ρf)');
ylabel('Water Amplitude (ρw)');
title('Optimization Success Rate (%)');
axis square;

saveas(gcf, 'fig/task7_stability_analysis.png');
