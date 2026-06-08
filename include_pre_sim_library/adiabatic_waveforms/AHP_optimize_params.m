function P = AHP_optimize_params(tau, wmax)

% Author: Maximilian Gram, University Hospital Wuerzburg, Wuerzburg, Germany; V1, 09.03.2026

    % define field properties for spins
    temp_db1 = 0.75 : 0.05 : 1.25;
    temp_df0 = -50  : 10   : 50;
    
    % define weights and convert to 1d
    temp_w1  = 1 - abs(linspace(-0.5, 0.5, numel(temp_db1)));
    temp_w2  = 1 - abs(linspace(-0.5, 0.5, numel(temp_df0)));
    c = 1;
    for j=1:numel(temp_df0)
    for k=1:numel(temp_db1)
        db1(c,1)    = temp_db1(k);
        df0(c,1)    = temp_df0(j);
        weight(c,1) = temp_w1(j) * temp_w2(k);
        c = c + 1;
    end
    end
    weight = weight / sum(weight);
    clear temp_db1 temp_df0 temp_w1 temp_w2 c j k;
    
    %% optimization   
    min_fun = @(P) get_cost(P, tau, wmax, db1, df0, weight);
    options = optimoptions('fmincon','Display','off','Algorithm','sqp');
    P       = fmincon(min_fun, [1, 1, 1], [], [], [], [], [0, 0, 0], [10, 10, 10], [], options);
    
    %% vis: compare initial vs optimized
    t = (0:1e-6:tau)';
    M = zeros(numel(db1), 3, numel(t));
    [w1, u] = AHP_modulation(t, tau, wmax, 1, 1, 1, 'tipdown');
    parfor j=1:numel(db1)
        M(j,:,:) = get_bloch_sim(w1, db1(j), df0(j));
    end
    figure()
    subplot(2,4,1)
    hold on
    plot(t*1e3, real(w1)/2/pi, 'LineWidth', 3)
    plot(t*1e3, imag(w1)/2/pi, 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['real/imag f1 [Hz]'])
    axis square
    title('initial')
    subplot(2,4,2)
    hold on
    plot(t*1e3, squeeze(M(:,1,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),1,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['Mx'])
    axis square
    subplot(2,4,3)
    hold on
    plot(t*1e3, squeeze(M(:,2,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),2,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['My'])
    axis square
    subplot(2,4,4)
    hold on
    plot(t*1e3, squeeze(M(:,3,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),3,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['Mz'])
    axis square
    
    t = (0:1e-6:tau)';
    M = zeros(numel(db1), 3, numel(t));
    [w1, u] = AHP_modulation(t, tau, wmax, P(1), P(2), P(3), 'tipdown');
    parfor j=1:numel(db1)
        M(j,:,:) = get_bloch_sim(w1, db1(j), df0(j));
    end
    subplot(2,4,5)
    hold on
    plot(t*1e3, real(w1)/2/pi, 'LineWidth', 3)
    plot(t*1e3, imag(w1)/2/pi, 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['real/imag f1 [Hz]'])
    axis square
    title('optimized')
    subplot(2,4,6)
    hold on
    plot(t*1e3, squeeze(M(:,1,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),1,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['Mx'])
    axis square
    subplot(2,4,7)
    hold on
    plot(t*1e3, squeeze(M(:,2,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),2,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['My'])
    axis square
    subplot(2,4,8)
    hold on
    plot(t*1e3, squeeze(M(:,3,:))')
    plot(t*1e3, squeeze(M(ceil(numel(db1)/2),3,:))', 'k-', 'LineWidth', 3)
    xlabel(['time [ms]'])
    ylabel(['Mz'])
    axis square

end

%% ------------ subroutines -------------
function [w1, u] = get_w1_u(P, t, tau, wmax)
    [w1, u] = AHP_modulation(t, tau, wmax, P(1), P(2), P(3), 'tipdown');
end

function M = get_bloch_sim(w1, db1, df0)
    n_dt = numel(w1);
    M    = zeros(3, n_dt);
    M_   = [0;0;1];
    for j=1:n_dt
        w1x_ = real(w1(j)) * db1;
        w1y_ = imag(w1(j)) * db1;
        dw0_ = 2*pi * df0;            
        B_   = [ 0,      dw0_,  -w1y_;
                -dw0_,   0,      w1x_;
                 w1y_,  -w1x_,   0];
        M_ = expm(B_*1e-6) * M_;
        M(:,j) = M_;
    end
end

function cost = get_cost(P, tau, wmax, db1, df0, weight)
    n       = 50;    
    t       = (0:1e-6:tau)';    
    n_dt    = numel(t);
    n_spins = numel(db1);
    cost    = zeros(n_spins, 1);
    [w1, u] = AHP_modulation(t, tau, wmax, P(1), P(2), P(3), 'tipdown');
    parfor j=1:n_spins
        M  = get_bloch_sim(w1, db1(j), df0(j));
        c1 = mean(sqrt(sum((M-u).^2))); % cost1: adiabatic criterion
        M  = M(:,end-n+1:end);        
        c2 = mean(sqrt((M(1,:)-1).^2 + M(2,:).^2 + M(3,:).^2)); % cost2: reaching [1;0;0]
        cost(j) = c1*c2;
    end
    cost = sum(cost.*weight);
end