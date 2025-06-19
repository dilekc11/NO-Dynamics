function NO_dynamics_simulation_finall()
    clear all; 
    
    params = struct(...
        'CaM_total', 1.0, ...    % uM (Total Calmodulin concentration)
        'K_d_CaM', 0.5, ...      % uM (Apparent dissociation constant for Ca2 to CaM)
        'n_Hill', 2.0, ...       % Hill coefficient for Ca2 binding to CaM
        'Vmax_nNOS', 25e-3, ... % uM/s 
        'Km_nNOS', 9.27e-2, ... % uM 
        'mu_deact_n', 0.0167, ... % s^-1 
        'Vmax_NO_n', 4.22, ...   % s^-1 
        'Km_O2_n', 243, ...      % uM 
        'KmLArg_n', 1.5, ...     % uM 
        'k_O2_n', 9.6e-6, ...    % uM^-2 s^-1 
        'X_nk', 25, ...          % um 
        'D_C_NO', 3300, ...      % um^2/s 
        'O2_n', 200, ...         % uM 
        'LArg_n', 100, ...       % uM 
        'V_spline', 8e-11, ...   % Dendritic spine volume (L)
        'K_ev', 1.6e3, ...       % s^-1 (Ca2 decay rate) 
        'Ca2_rest', 1e-7, ...   % uM (resting [Ca2])
        'A_buf', 20, ...         % Buffer capacity (dimensionless) 
        'v_n', -0.04, ...       % Neuronal membrane potential (V) 
        'G_M', 4.6e4, ...       % NMDA-R conductance (fS). 
        'P_Ca_over_P_M', 3.6, ...% Permeability ratio (P_Ca/P_M). 
        'Ca_ex', 2000, ...       % External [Ca2] (uM) 
        'M', 1.3e5, ...         % Monovalent ion concentration (uM) 
        'alpha_v', -80, ...     % Voltage translation factor (V^-1) 
        'beta_v', 0.02, ...      % Voltage translation offset (V) 
        'T', 310, ...            % Temperature (K) 
        'R_gas', 8.314, ...      % Gas constant (J/(mol.K)) 
        'F', 96485, ...          % Faraday constant (C/mol) 
        'Km_A', 650, ...         % Michaelis constant for NR2A (uM) 
        'Km_B', 2800, ...         % Michaelis constant for NR2B (uM) 
        't0', 400, ...            % Start time for Glu pulse (s) 
        't2', 600, ...            % End time for Glu pulse (s)  
        'Glu_max', 1846 ...     % Maximum glutamate concentration (uM)
    );
    
    fprintf("\n--- Debugging Information ---\n");
    fprintf("  Ca2_n_initial (set to params.Ca2_rest): %.7f uM\n", params.Ca2_rest); 
    fprintf("  params.Ca2_rest: %.7f uM\n", params.Ca2_rest); 
    fprintf("  params.K_d_CaM: %.2f uM\n", params.K_d_CaM);
    fprintf("  params.n_Hill: %.2f\n", params.n_Hill);
    fprintf("  params.CaM_total: %.2f uM\n", params.CaM_total);
    fprintf("-----------------------------\n\n");

    Ca2_n_initial = params.Ca2_rest; 
    tspan = [0, 1000];         
    y0 = [Ca2_n_initial, 0, 0, 0]; 
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-13);

    [t, y] = ode45(@(t, y) no_dynamics(t, y, params), tspan, y0, options);

    function dydt = no_dynamics(t, y, params)
        Ca2_n = y(1);
        nNOS_scL = y(2);
        NO_n = y(3);
        NO_k = y(4);

        %% 1. Glutamate
        if t >= params.t0 && t <= params.t2
            Glul_sc = params.Glu_max * (0.5 * tanh(t - params.t0) - 0.5 * tanh(t - params.t2));
        else
            Glul_sc = 0;
        end
        
        %% 2. NMDA Receptor Activation and inward calcim current
        w_NR2A = Glul_sc / (params.Km_A + Glul_sc);
        w_NR2B = Glul_sc / (params.Km_B + Glul_sc);

        numerator = 4 * params.v_n * params.G_M * params.P_Ca_over_P_M * (params.Ca_ex / params.M);
        denominator = 1 + exp(params.alpha_v * (params.v_n + params.beta_v));
        exponent_term = exp(2 * params.v_n * params.F / (params.R_gas * params.T));
        fraction = exponent_term / (1 - exponent_term);
        
        I_Ca_fA = (numerator / denominator) * fraction;
        I_Ca_tot_fA = I_Ca_fA * (0.63 * w_NR2A + 11 * w_NR2B);
        I_Ca_tot_A = I_Ca_tot_fA * 1e-15;
        
        %% 3. Influx Ca2+
        term1 = (-I_Ca_tot_A / (2 * params.F * params.V_spline)) * 1e6;
        term2 = params.K_ev * (Ca2_n - params.Ca2_rest);
        dCa2_n_dt = (term1 - term2) / (1 + params.A_buf);

        %% 4. Ca-CaM Complex Contentration
        CaM_n = params.CaM_total * (Ca2_n^params.n_Hill) / (params.K_d_CaM^params.n_Hill + Ca2_n^params.n_Hill);

        %% 5. nNOS Activation
        dnNOS_scL_dt = (params.Vmax_nNOS * CaM_n) / (params.Km_nNOS + CaM_n) - params.mu_deact_n * nNOS_scL;

        %% 6. NO Proudction
        p_NO_n = params.Vmax_NO_n * nNOS_scL * ...
                 (params.O2_n / (params.Km_O2_n + params.O2_n)) * ...
                 (params.LArg_n / (params.KmLArg_n + params.LArg_n));

        %7 NO Consumption Flux
        c_NO_n = params.k_O2_n * NO_n^2 * params.O2_n;
        %Time for NO to diffuse between the centres of the NE and the AC
        tau_ns = (params.X_nk^2) / (2 * params.D_C_NO);
        %NO Diffusion Flux
        d_NO_n = (NO_k - NO_n) / tau_ns;
        %8  
        dNO_n_dt = p_NO_n - c_NO_n + d_NO_n; 

        dydt = [dCa2_n_dt; dnNOS_scL_dt; p_NO_n - c_NO_n + d_NO_n; 0]; 
    end

    Glu_conc = params.Glu_max * (0.5 * tanh(t - params.t0) - 0.5 * tanh(t - params.t2));
    w_NR2A_plot = Glu_conc ./ (params.Km_A + Glu_conc);
    w_NR2B_plot = Glu_conc ./ (params.Km_B + Glu_conc);

    % Calculate CaM_n for plotting (using the same formula as in the ODE)
    CaM_n_plot = params.CaM_total * (y(:,1).^params.n_Hill) ./ (params.K_d_CaM^params.n_Hill + y(:,1).^params.n_Hill);

    figure;
    
    %% Plots
    subplot(7,1,1);
    plot(t, Glu_conc, 'LineWidth', 2, 'Color', 'g');
    xlabel('Time (s)');
    ylabel('[Glu] (uM)');
    title('Glutamate Concentration at synaptic cleft');
    grid on;

    subplot(7,1,2);
    plot(t, w_NR2A_plot, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    xlabel('Time (s)');
    ylabel('w_NR2A');
    title('NR2A Weight (NMDA Receptor)');
    grid on;

    subplot(7,1,3);
    plot(t, w_NR2B_plot, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
    xlabel('Time (s)');
    ylabel('w_NR2B');
    title('NR2B Weight (NMDA Receptor)');
    grid on;

    subplot(7,1,4);
    plot(t, y(:,1), 'LineWidth', 2, 'Color', 'r');
    xlabel('Time (s)');
    ylabel('[Ca2+]_n (uM)');
    title('Calcium Concentration Dynamics (Influx Ca2+)');
    grid on;

    subplot(7,1,5);
    plot(t, CaM_n_plot, 'LineWidth', 2, 'Color', 'm');
    xlabel('Time (s)');
    ylabel('[CaM]_n (uM)');
    title('CaM Concentration Dynamics');
    grid on;

    subplot(7,1,6);
    plot(t, y(:,2), 'LineWidth', 2, 'Color', 'b');
    xlabel('Time (s)');
    ylabel('[nNOS_act] (uM)');
    title('nNOS Activation');
    grid on;

    subplot(7,1,7);
    plot(t, y(:,3), 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('[NO]_n (uM)');
    title('NO Concentration Dynamics');
    grid on;
    
    figure;
    plot(y(:,1), y(:,3), 'LineWidth', 2);
    xlabel('[Ca2+]_n (uM)');
    ylabel('[NO]_n (uM)');
    title('NO Concentration vs Calcium Concentration');
    grid on;
    
    steady_state_NO = mean(y(end-10:end, 3));
    steady_state_Ca = mean(y(end-10:end, 1));
    fprintf('Steady-state values:\n');
    fprintf('  [NO]_n = %.4f uM\n', steady_state_NO);
    fprintf('  [Ca2+]_n = %.7f uM\n', steady_state_Ca); 
end


