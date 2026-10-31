%Task 1

T    = (0:0.001:0.010)';
am_w = 1;          % water amplitude
am_f = 0.5;        % fat amplitude
R2s   = 50;         % R2* decay (s^-1)
ff    = -0.45e3;    % single‐peak fat freq (Hz)
r    = [0.087,0.694,0.128,0.004,0.039,0.048]';
fF_m = [-0.498,-0.447,-0.345,-0.261,-0.063,0.064]*1e3;

%--- Mi: single-peak, no decay
S1 = am_w + am_f * exp(1i*2*pi*ff*T);

%--- Mii: multi-peak, no decay
S2 = am_w + am_f * ( exp(1i*2*pi*(T * fF_m)) * r );
%--- Miii: multi-peak + R2* decay
S3 = S2 .* exp(-T * R2s);

% Start plotting
fig = figure('Name','Signal Evolution Models','NumberTitle','off','Position',[100 100 900 600]);

% 1) Mag vs Time
subplot(2,2,1);
plot(T*1e3, abs(S1),'-','LineWidth',1.5,'Color',[0 0.4470 0.7410]); hold on;
plot(T*1e3, abs(S2),'--','LineWidth',1.5,'Color',[0.8500 0.3250 0.0980]);
plot(T*1e3, abs(S3),':','LineWidth',1.5,'Color',[0.4660 0.6740 0.1880]);
hold off;
xlabel('Time (ms)');
ylabel('Signal Magnitude');
legend('Model (i)','Model (ii)','Model (iii)','Location','northeast');
title('Signal Magnitude over Time');
grid on;

% 2) Real  vs Time
subplot(2,2,2);
plot(T*1e3, real(S1),'-o','MarkerSize',4); hold on;
plot(T*1e3, real(S2),'-s','MarkerSize',4);
plot(T*1e3, real(S3),'-^','MarkerSize',4);
hold off;
xlabel('Time (ms)');
ylabel('Real Part');
legend('Model (i)','Model (ii)','Model (iii)','Location','northeast');
title('Real Part over Time');
grid on;

% 3) Imag  vs Time
subplot(2,2,3);
plot(T*1e3, imag(S1),'-o','MarkerSize',4); hold on;
plot(T*1e3, imag(S2),'-s','MarkerSize',4);
plot(T*1e3, imag(S3),'-^','MarkerSize',4);
hold off;
xlabel('Time (ms)');
ylabel('Imag Part');
legend('Model (i)','Model (ii)','Model (iii)','Location','northeast');
title('Imag Part over Time');
grid on;

% 4) Complex trajec
subplot(2,2,4);
plot(real(S1), imag(S1),'-o','MarkerSize',4); hold on;
plot(real(S2), imag(S2),'-s','MarkerSize',4);
plot(real(S3), imag(S3),'-^','MarkerSize',4);
hold off;
xlabel('Real Part');
ylabel('Imag Part');
legend('Model (i)','Model (ii)','Model (iii)','Location','northeast');
title('Signal in Complex Plane');
axis equal;
grid on;

sgtitle('Signal Evolution Models: Magnitude, Real, Imag & Complex Plane');
saveas(fig, 'fig/task1.png');



% Task 2
sets_two= [1,0.1;  0.1,1;  1,0.3;  0.3,1];
fatSignal = exp(1i*2*pi*(T*fF_m)) * r;

fig = figure('Name','Task 2: Fat–Water Ambiguity','NumberTitle','off','Position',[200 200 800 600]);

for k = 1:4
    rw = sets_two(k,1);
    rf = sets_two(k,2);

    % Model (i)
    S_i   = rw + rf * exp(1i*2*pi*ff*T);
    % Model (ii)
    S_ii  = rw + rf * fatSignal ;
    % Model (iii)
    S_iii = S_ii .* exp(-T*R2s);

    % Plot magnitudes
    subplot(2,2,k);
    plot(T*1e3, abs(S_i),  '-o','LineWidth',1.2,'MarkerSize',4,'DisplayName','(i)'); hold on;
    plot(T*1e3, abs(S_ii), '--s','LineWidth',1.2,'MarkerSize',4,'DisplayName','(ii)');
    plot(T*1e3, abs(S_iii),':^','LineWidth',1.2,'MarkerSize',4,'DisplayName','(iii)');
    hold off;

    xlabel('T (ms)');
    ylabel('|S|');
    title(sprintf('am_w=%.1f, am_f=%.1f', rw, rf));
    legend('Location','northeast');
    grid on;
end

sgtitle('Task 2: Magnitude vs T for Different am_w,am_f under Models (i),(ii),(iii)');
saveas(fig, 'fig/task2.png');


% Task 3
SNR = 70;
% from Task 2
sets_two = [1.0, 0.1;
             0.1, 1.0;
             1.0, 0.3;
             0.3, 1.0];
R2s_true = 70;


Nstarts = 20;
Nwater = ceil(Nstarts/2);
Nfat = Nstarts - Nwater;
lb = [0, 0, 50];
ub = [1, 1, 90];

opts = optimoptions('fmincon','Display','off');
fig3 = figure('Name','Parameter Estimation Analysis','NumberTitle','off','Position',[100 100 1200 800]);

for test_idx = 1:4
    am_w_true = sets_two(test_idx,1);
    am_f_true = sets_two(test_idx,2);

    fatSignal = exp(1i*2*pi*(T*fF_m)) * r;
    S_true = (am_w_true + am_f_true*fatSignal) .* exp(-T*R2s_true);
    noisyMag = abs(S_true + (max(abs(S_true))/SNR)*(randn(size(S_true))+1i*randn(size(S_true))));

    all_params = zeros(Nstarts, 3);
    all_SSE = zeros(Nstarts, 1);
    init_params = zeros(Nstarts, 3);
    % sse
    costfun = @(p) sum((noisyMag - abs((p(1) + p(2)*fatSignal) .* exp(-T*p(3)))).^2);


    for k = 1:Nwater
        if k == 1
            init = [1, 0, 70];
        elseif k == 2
            init = [0.9, 0.1, 65];
        elseif k == 3
            init = [0.7, 0.3, 75];
        else
            init = [0.7 + 0.2*rand, 0.1 + 0.2*rand, 65 + 10*rand];
        end
        init_params(k,:) = init;

        [p,ss] = fmincon(costfun, init, [],[],[],[], lb, ub, [], opts);
        all_params(k,:) = p;
        all_SSE(k) = ss;
    end

    for k = 1:Nfat
        if k == 1
            init = [0, 1, 70];
        elseif k == 2
            init = [0.1, 0.9, 65];
        elseif k == 3
            init = [0.3, 0.7, 75];
        else
            init = [0.1 + 0.2*rand, 0.7 + 0.2*rand, 65 + 10*rand];
        end
        init_params(Nwater+k,:) = init;

        [p,ss] = fmincon(costfun, init, [],[],[],[], lb, ub, [], opts);
        all_params(Nwater+k,:) = p;
        all_SSE(Nwater+k) = ss;
    end

    [best_SSE, best_idx] = min(all_SSE);
    best_params = all_params(best_idx,:);

    % cluster analysis (guard against division by zero for small params)
    tolerance = 0.05;
    denom = max(abs(best_params), 1e-6);
    similar_to_best = sum(abs(all_params - best_params) ./ denom < tolerance, 2) == 3;
    success_rate = sum(similar_to_best) / Nstarts * 100;

    subplot(2,2,test_idx);
    plot(T*1e3, noisyMag, 'ko', 'MarkerSize', 4, 'DisplayName', 'Noisy |S|'); hold on;
    S_best = (best_params(1) + best_params(2)*fatSignal) .* exp(-T*best_params(3));
    plot(T*1e3, abs(S_best), 'r-', 'LineWidth', 1.5, 'DisplayName', 'Best Fit');
    plot(T*1e3, abs(S_true), 'b--', 'LineWidth', 1.5, 'DisplayName', 'True Signal');
    xlabel('T (ms)'); ylabel('|S|');
    legend('Location', 'northeast');
    title(sprintf('Test Set %d (ρw=%.1f, ρf=%.1f)\nSSE=%.4f', test_idx, am_w_true, am_f_true, best_SSE));
    grid on;

    fprintf('\nTest Set %d Analysis (ρw=%.1f, ρf=%.1f):\n', test_idx, am_w_true, am_f_true);
    fprintf('Success rate (within 5%% of best): %.1f%%\n', success_rate);
    fprintf('Best fit parameters (SSE = %.4f):\n', best_SSE);
    fprintf('  Water amp: %.4f (true: %.4f, error: %.1f%%)\n', best_params(1), am_w_true, abs(best_params(1)-am_w_true)/am_w_true*100);
    fprintf('  Fat amp:  %.4f (true: %.4f, error: %.1f%%)\n', best_params(2), am_f_true, abs(best_params(2)-am_f_true)/am_f_true*100);
    fprintf('  R2*:      %.2f (true: %.2f, error: %.1f%%)\n', best_params(3), R2s_true, abs(best_params(3)-R2s_true)/R2s_true*100);
end

sgtitle('Task 3: Model Fitting Results for All Test Sets');
saveas(fig3, 'fig/task3_analysis.png');
