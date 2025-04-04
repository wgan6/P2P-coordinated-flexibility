clear;clc;
close all;
T = 24;  % Number of time slots (e.g., hours in a day)

% Regional scaling factors for PV penetration rates, according to FES
regionRate = [1.4 1.6 2.4 3.1 1.4 1.2 1.2 1.2];

% PV penetration base values for LW and CT scenarios over time (2020–2050)
PVpene_rateSet0 = [0.05 0.08 0.15 0.23 0.32 0.39 0.44 0.05 0.08 0.15 0.23 0.32 0.39 0.44];

% ESS penetration rate for each year/scenario (CT capacity is 87% of LW)
ESSrateSet = [0.00 0.01 0.05 0.08 0.11 0.14 0.18 0.00 0.01 0.05 0.08 0.11 0.14 0.18];

% DR (Demand Response) penetration rates
DR_rateSet = [0.00 0.00 0.08 0.18 0.23 0.23 0.23 0.00 0.00 0.06 0.16 0.21 0.21 0.21];

% Optimization tolerance and pricing switch
Permit_costDeviation = 1e-3;  % Allowable deviation in cost objectives
ToUPri = 1;                         % Time-of-Use pricing flag


regionInd = [1]; % selected region

% Adjust PV penetration values according to regional scaling factor
PVpene_rateSet = roundn(PVpene_rateSet0 * regionRate(regionInd) / regionRate(1), -2);

% Select scenario index (e.g., 7 corresponds to LW2050)
rateIIInd = [7];
ESSrateInd = [rateIIInd];
PVpene_rateInd = [rateIIInd];
DR_rateInd = [rateIIInd];

% Choose specific load and PV scenario
LDsceInd = [13];
PVsceInd = [3];

% Select network configuration
Network_ID = [1];

% Define scenario variation set
IIIset2 = [2001];
IIIset1 = IIIset2 - 1000;
IIIset3 = IIIset2 + 1000;
IIIset123 = reshape([IIIset2; IIIset1; IIIset3], 1, length(IIIset1) * 3); % P2P, No P2P, No Flex,
IIIset12 = reshape([IIIset2; IIIset1], 1, length(IIIset1) * 2);

for III = IIIset123

    % Load network, load, PV, EV and other input data for the scenario
    Case_Parameters_from_UK_FES;

    tic;  % Start timing the optimization

    %% ---------- Variable Declarations ----------
    %% Power system variables
    Ps = sdpvar(Ns, T, 'full');         % Substation active power
    Qs = sdpvar(Ns, T, 'full');         % Substation reactive power
    Pp2p_sell = sdpvar(1, T, 'full');   % Power sold in P2P market
    Pp2p_buy = sdpvar(1, T, 'full');    % Power bought from P2P market
    Pl = sdpvar(Nl, T, 'full');         % Line active power
    Ql = sdpvar(Nl, T, 'full');         % Line reactive power
    Usqa = sdpvar(Nb, T, 'full');       % Squared voltage magnitude at buses
    PV_out = sdpvar(Nb, T, 'full');     % PV generation output

    %% ESS (Energy Storage System) variables
    Pes_ch = sdpvar(Nb, T, 'full');     % Charging power of ESS
    Pes_dis = sdpvar(Nb, T, 'full');    % Discharging power of ESS
    Pes = sdpvar(Nb, T, 'full');        % Net ESS power
    Ees = sdpvar(Nb, T, 'full');        % Energy level of ESS
    Eini = sdpvar(Nb, 1, 'full');       % Initial energy level of ESS

    %% V2G (Vehicle-to-Grid) energy storage variables
    Pes_ch_V2G = sdpvar(EVnum_V2G, T, 'full');  % V2G charging power
    Pes_dis_V2G = sdpvar(EVnum_V2G, T, 'full'); % V2G discharging power
    Pes_V2G = sdpvar(EVnum_V2G, T, 'full');     % Net V2G power
    Ees_V2G = sdpvar(EVnum_V2G, T, 'full');     % V2G energy level

    %% Demand Response variables
    delpi0 = sdpvar(Nb, T, 'full');     % Change in baseline demand
    pi0Actu = sdpvar(Nb, T, 'full');    % Actual load after DR
    PnetLD = sdpvar(Nb, T, 'full');     % Net load demand
    Pfeed = sdpvar(Nb, T, 'full');      % Power fed into the grid

    %% ---------- Power System Constraints ----------
    Constraints = [];

    % Demand Response (DR) constraints: bounds and conservation
    Constraints = [Constraints; pi0min <= pi0Actu <= pi0max];         % DR limits per node and time
    Constraints = [Constraints; sum(pi0Actu, 2) == sum(pi0, 2)];      % Energy conservation (same daily total)
    Constraints = [Constraints; delpi0 >= pi0Actu - pi0];            % DR shift must match difference
    Constraints = [Constraints; delpi0 >= 0];                        % DR must be non-negative

    % ---------------------- Power Flow & Balance Constraints ----------------------
    % Ensure node-level power balance: generation, PV, ESS, DR, and EVs equal demand
    Constraints = [Constraints; Aim' * Pl == [Ps; zeros(Nb - 1, T)] + PV_out - Pes - pi0Actu - Pev_Dumb_Smart - Hev * Pes_V2G];
    Constraints = [Constraints; Aim' * Ql == [Qs; zeros(Nb - 1, T)] - qi0];
    Constraints = [Constraints; PnetLD - Pfeed == pi0Actu + Pev_Dumb_Smart + Hev * Pes_V2G + Pes - PV_out];

    % PV generation and power boundaries
    Constraints = [Constraints; 0 <= PV_out <= PV_pre];  % PV output within forecasted availability
    Constraints = [Constraints; 0 <= PnetLD];           % Net demand must be positive
    Constraints = [Constraints; 0 <= Pfeed];            % Feed-in power to grid must be positive

    % ---------------------- ESS Constraints ----------------------
    % Energy Storage System power and energy flow constraints
    Constraints = [Constraints; Pes == Pes_ch - Pes_dis];                                % Net power flow
    Constraints = [Constraints; 0 <= Pes_ch <= repmat(cmax, 1, T)];                      % Charging power bounds
    Constraints = [Constraints; 0 <= Pes_dis <= repmat(cmax, 1, T)];                     % Discharging power bounds
    Constraints = [Constraints; 0 <= Ees <= repmat(smax, 1, T)];                         % Energy bounds
    Constraints = [Constraints; Ees(:, T) == Eini];                                      % ESS must return to initial energy at end of horizon
    Constraints = [Constraints; Ees == [Eini, Ees(:, 1:T - 1)] + Pes_ch * eff_es - Pes_dis / eff_es];  % Energy transition

    % ---------------------- V2G Constraints (per Electric Vehicle) ----------------------
    % Rules for when EVs are connected or disconnected, and their energy states
    for Ind = 1:EVnum_V2G
        if ismember(Ind, evNoTranver)
            Constraints = [Constraints; Pes_V2G(Ind, EV_outTime(Ind)) == 0];
        elseif ismember(Ind, evTranver)
            Constraints = [Constraints; Pes_V2G(Ind, EV_adTime(Ind).one(end)) == 0];
        end
        Constraints = [Constraints; Pes_V2G(Ind, EV_disconnecTime(Ind).one) == 0];
        Constraints = [Constraints; Ees_V2G(Ind, EV_disconnecTime(Ind).one) == 0];
        Constraints = [Constraints; 0 <= Ees_V2G(Ind, EV_connecTime(Ind).one) <= repmat(smax_V2G(Ind), 1, length(EV_connecTime(Ind).one))];
    end

    % V2G net power = charging - discharging
    Constraints = [Constraints; Pes_V2G == Pes_ch_V2G - Pes_dis_V2G];
    Constraints = [Constraints; 0 <= Pes_ch_V2G <= repmat(cmax_V2G, 1, T)];
    Constraints = [Constraints; 0 <= Pes_dis_V2G <= repmat(cmax_V2G, 1, T)];

    % V2G energy transition constraints (non-travel and travel cases)
    for Ind = 1:EVnum_V2G
        Tin = EV_inTime(Ind);      % Start charging
        Tout = EV_outTime(Ind);    % End charging (within 24h)

        if ismember(Ind, evNoTranver)
            Constraints = [Constraints; Ees_V2G(Ind, Tin) == Ein_V2G(Ind)];
            Constraints = [Constraints; Ees_V2G(Ind, Tin + 1:Tout) == Ees_V2G(Ind, Tin:Tout - 1) + Pes_ch_V2G(Ind, Tin:Tout - 1) * eff_es - Pes_dis_V2G(Ind, Tin:Tout - 1) / eff_es];
            Constraints = [Constraints; Ees_V2G(Ind, Tout) == Eout_V2G(Ind)];
        end

        if ismember(Ind, evTranver)
            Constraints = [Constraints; Ees_V2G(Ind, Tin) == Ein_V2G(Ind)];
            Constraints = [Constraints; Ees_V2G(Ind, Tin + 1:24) == Ees_V2G(Ind, Tin:23) + Pes_ch_V2G(Ind, Tin:23) * eff_es - Pes_dis_V2G(Ind, Tin:23) / eff_es];
            Constraints = [Constraints; Ees_V2G(Ind, 1) == Ees_V2G(Ind, 24) + Pes_ch_V2G(Ind, 24) * eff_es - Pes_dis_V2G(Ind, 24) / eff_es];

            Tadd = EV_adTime(Ind).one;
            if Tadd(end) >= 2
                Constraints = [Constraints; Ees_V2G(Ind, 2:Tadd(end)) == Ees_V2G(Ind, 1:Tadd(end) - 1) + Pes_ch_V2G(Ind, 1:Tadd(end) - 1) * eff_es - Pes_dis_V2G(Ind, 1:Tadd(end) - 1) / eff_es];
            end
            Constraints = [Constraints; Ees_V2G(Ind, Tadd(end)) == Eout_V2G(Ind)];
        end
    end

    % ---------------------- Voltage and System Limit Constraints ----------------------
    % Linearized voltage drop model across the distribution network
    Constraints = [Constraints; Aim * Usqa == 2 * Pl .* repmat(Dis.Rexi, 1, T) + 2 * Ql .* repmat(Dis.Xexi, 1, T)];

    % Physical limits on line power, substation output, and voltages
    for t = 1:T
        for i = 1:Nl
            Constraints = [Constraints; cone([Dis.Ecap(i); Pl(i,t); Ql(i,t)])];
        end
    end
    for t = 1:T
        for i = 1:Ns
            Constraints = [Constraints; cone([Scap; Ps(i,t); Qs(i,t)])];
        end
    end
    Constraints = [Constraints; Vi_min^2 <= Usqa <= Vi_max^2];

    % ---------------------- P2P Market Energy Balance ----------------------
    % Ensure peer-to-peer market has equal buy/sell balance over the network
    Constraints = [Constraints; Pp2p_buy >= 0];
    Constraints = [Constraints; Pp2p_sell >= 0];
    Constraints = [Constraints; sum(PnetLD(p2pSet, :) - Pfeed(p2pSet, :), 1) == Pp2p_buy - Pp2p_sell];

    %% ---------- First-Stage Objective Function and Optimization ----------

    % Cost of Demand Response (DR) = total shifted load × unit DR cost
    Fdr = sum(sum(delpi0)) * DR_cost;

    % Cost of ESS degradation = (charging + discharging energy) × unit ESS cost × degradation factor
    Fes = sum(sum(k_degra / 100 * (Pes_dis + Pes_ch) .* ESS_capcost_perunit .* repmat(EShour_duration, 1, T)));

    % Cost of V2G ESS degradation (same as above but for V2G batteries)
    Fes_V2G = sum(sum(k_degra / 100 * (Pes_dis_V2G + Pes_ch_V2G) .* ESS_capcost_perunit_V2G .* repmat(EShour_duration_V2G, 1, T)));

    % Cost of electricity consumption from grid: buying at ToU price minus selling via FiT
    Fpdn1 = sum(rou_ToU .* sum(PnetLD(singleSet, :), 1)) - sum(rou_FiT .* sum(Pfeed(singleSet, :), 1));

    % Cost/revenue from P2P trading
    Fpdn2 = sum(rou_ToU .* Pp2p_buy) - sum(rou_FiT .* Pp2p_sell);

    % Total cost as the optimization objective (to be minimized)
    Objective = Fpdn1 + Fpdn2 + Fes + Fes_V2G + Fdr;

    % Solver settings using Gurobi, with verbose output
    options = sdpsettings('solve', 'mosek');
    options = sdpsettings(options, 'verbose', 2);

    % Solve the first optimization problem
    solvesdp(Constraints, Objective, options);

    % ---------- Save Initial Solution Results for Further Use ----------

    Fpdn_init = double(Objective);
    Fpdn1_init = double(Fpdn1);
    Fpdn2_init = double(Fpdn2);
    Fes_init = double(Fes);
    Fes_V2G_init = double(Fes_V2G);
    Fdr_init = double(Fdr);

    % P2P hourly trading cost/revenue
    Fpdn2Hourly_init = double(rou_ToU .* Pp2p_buy - rou_FiT .* Pp2p_sell);

    % Store optimized decision variable values for second-stage reuse
    Pfeed_init = double(Pfeed);
    PnetLD_init = double(PnetLD);

    pi0Actu_init = double(pi0Actu);
    Pes_init = double(Pes);
    PV_out_init = double(PV_out);

    Pp2p_buy_init = double(Pp2p_buy);
    Pp2p_sell_init = double(Pp2p_sell);
    Ees_init = double(Ees);

    Pes_ch_V2G_init = double(Pes_ch_V2G);
    Pes_dis_V2G_init = double(Pes_dis_V2G);
    Pes_V2G_init = double(Pes_V2G);
    Ees_V2G_init = double(Ees_V2G);

    Ps_init = double(Ps);
    Pl_init = double(Pl);
    Ql_init = double(Ql);

    % ---------- Asset Life Expectancy-Based Investment Evaluation ----------

    % Calculate Substation and Line Spare Capacity Ratio (SCR)
    SCRL_init = (Dis.Ecap - max(abs(Pl_init), [], 2)) ./ Dis.Ecap;
    SCRS_init = (Scap - max(abs(Ps_init), [], 2)) ./ Scap;

    % Estimate Present Value of Line (PVL) and Substation (PVS) assets
    PVL_init = AssetL .* (1 + discountR).^(log(1 - SCRL_init) / log(1 + loadR));
    PVS_init = AssetS * (1 + discountR).^(log(1 - SCRS_init) / log(1 + loadR));

    % Aggregate total asset present value
    PVA_init = PVS_init + sum(PVL_init);


    %%                Further Optimization
    %%                Further Optimization
    %% ----------------------------------------------------
    % Start second-stage optimization if P2P trading is enabled
    if p2pRate > 0

        % --------- Redefine Variables for Full Re-optimization ---------
        % Power system variables (same structure as in stage 1)
        Ps = sdpvar(Ns, T, 'full');
        Qs = sdpvar(Ns, T, 'full');
        Pp2p_sell = sdpvar(1, T, 'full');
        Pp2p_buy = sdpvar(1, T, 'full');
        Pl = sdpvar(Nl, T, 'full');
        Ql = sdpvar(Nl, T, 'full');
        Usqa = sdpvar(Nb, T, 'full');
        PV_out = sdpvar(Nb, T, 'full');

        % ESS variables
        Pes_ch = sdpvar(Nb, T, 'full');
        Pes_dis = sdpvar(Nb, T, 'full');
        Pes = sdpvar(Nb, T, 'full');
        Ees = sdpvar(Nb, T, 'full');
        Eini = sdpvar(Nb, 1, 'full');

        % V2G variables
        Pes_ch_V2G = sdpvar(EVnum_V2G, T, 'full');
        Pes_dis_V2G = sdpvar(EVnum_V2G, T, 'full');
        Pes_V2G = sdpvar(EVnum_V2G, T, 'full');
        Ees_V2G = sdpvar(EVnum_V2G, T, 'full');

        % Demand response variables
        delpi0 = sdpvar(Nb, T, 'full');
        pi0Actu = sdpvar(Nb, T, 'full');
        PnetLD = sdpvar(Nb, T, 'full');
        Pfeed = sdpvar(Nb, T, 'full');

        % ---------- SCR & Asset Present Value Variables ----------
        % Introduce virtual capacity variables and related marginal cost terms
        VirCapL = sdpvar(Nl, 1, 'full');        % Virtual line capacity
        SCRL = sdpvar(Nl, 1, 'full');           % Line spare capacity ratio
        delMarL = sdpvar(Nl, step, 'full');     % Piecewise margin adjustments for lines
        PVL = sdpvar(Nl, 1, 'full');            % Present value of line assets

        VirCapS = sdpvar(Ns, 1, 'full');        % Virtual substation capacity
        SCRS = sdpvar(Ns, 1, 'full');           % Substation spare capacity ratio
        delMarS = sdpvar(Ns, step, 'full');     % Piecewise margin adjustments for substations
        PVS = sdpvar(Ns, 1, 'full');            % Present value of substation assets

        %% ---------- Power System Constraints ----------
        Constraints = [];

        % Demand Response constraints (same as in stage 1)
        Constraints = [Constraints; pi0min <= pi0Actu <= pi0max];
        Constraints = [Constraints; sum(pi0Actu, 2) == sum(pi0, 2)];
        Constraints = [Constraints; delpi0 >= pi0Actu - pi0];
        Constraints = [Constraints; delpi0 >= 0];

        % Power balance equations (active/reactive)
        Constraints = [Constraints; Aim' * Pl == [Ps; zeros(Nb - 1, T)] + PV_out - Pes - pi0Actu - Pev_Dumb_Smart - Hev * Pes_V2G];
        Constraints = [Constraints; Aim' * Ql == [Qs; zeros(Nb - 1, T)] - qi0];
        Constraints = [Constraints; PnetLD - Pfeed == pi0Actu + Pev_Dumb_Smart + Hev * Pes_V2G + Pes - PV_out];

        % PV and net demand limits
        Constraints = [Constraints; 0 <= PV_out <= PV_pre];
        Constraints = [Constraints; 0 <= PnetLD];
        Constraints = [Constraints; 0 <= Pfeed];

        % ESS operation constraints (charging/discharging/energy)
        Constraints = [Constraints; Pes == Pes_ch - Pes_dis];
        Constraints = [Constraints; 0 <= Pes_ch <= repmat(cmax, 1, T)];
        Constraints = [Constraints; 0 <= Pes_dis <= repmat(cmax, 1, T)];
        Constraints = [Constraints; 0 <= Ees <= repmat(smax, 1, T)];
        Constraints = [Constraints; Ees == [Eini, Ees(:, 1:T-1)] + Pes_ch * eff_es - Pes_dis / eff_es];
        Constraints = [Constraints; Ees(:, T) == Eini];  % Final energy level same as initial

        % ---------- V2G (Vehicle-to-Grid) Constraints ----------
        for Ind = 1:EVnum_V2G
            % Ensure no discharging at disconnection time
            if ismember(Ind, evNoTranver)
                Constraints = [Constraints; Pes_V2G(Ind, EV_outTime(Ind)) == 0];
            elseif ismember(Ind, evTranver)
                Constraints = [Constraints; Pes_V2G(Ind, EV_adTime(Ind).one(end)) == 0];
            end
            % Force zero output and a dummy -10 SOC at disconnect time
            Constraints = [Constraints; Pes_V2G(Ind, EV_disconnecTime(Ind).one) == 0];
            Constraints = [Constraints; Ees_V2G(Ind, EV_disconnecTime(Ind).one) == -10];
            % Ensure SOC within limits during connection
            Constraints = [Constraints; 0 <= Ees_V2G(Ind, EV_connecTime(Ind).one) <= repmat(smax_V2G(Ind), 1, length(EV_connecTime(Ind).one))];
        end

        % Define net V2G power and bounds
        Constraints = [Constraints; Pes_V2G == Pes_ch_V2G - Pes_dis_V2G];
        Constraints = [Constraints; 0 <= Pes_ch_V2G <= repmat(cmax_V2G, 1, T)];
        Constraints = [Constraints; 0 <= Pes_dis_V2G <= repmat(cmax_V2G, 1, T)];

        % V2G energy transition modeling
        for Ind = 1:EVnum_V2G
            Tin = EV_inTime(Ind);
            Tout = EV_outTime(Ind);

            if ismember(Ind, evNoTranver)
                Constraints = [Constraints; Ees_V2G(Ind, Tin) == Ein_V2G(Ind)];
                Constraints = [Constraints; Ees_V2G(Ind, Tin+1:Tout) == Ees_V2G(Ind, Tin:Tout-1) + Pes_ch_V2G(Ind, Tin:Tout-1)*eff_es - Pes_dis_V2G(Ind, Tin:Tout-1)/eff_es];
                Constraints = [Constraints; Ees_V2G(Ind, Tout) == Eout_V2G(Ind)];
            end

            if ismember(Ind, evTranver)
                Constraints = [Constraints; Ees_V2G(Ind, Tin) == Ein_V2G(Ind)];
                Constraints = [Constraints; Ees_V2G(Ind, Tin+1:24) == Ees_V2G(Ind, Tin:23) + Pes_ch_V2G(Ind, Tin:23)*eff_es - Pes_dis_V2G(Ind, Tin:23)/eff_es];
                Constraints = [Constraints; Ees_V2G(Ind, 1) == Ees_V2G(Ind, 24) + Pes_ch_V2G(Ind, 24)*eff_es - Pes_dis_V2G(Ind, 24)/eff_es];

                Tadd = EV_adTime(Ind).one;
                if Tadd(end) >= 2
                    Constraints = [Constraints; Ees_V2G(Ind, 2:Tadd(end)) == Ees_V2G(Ind, 1:Tadd(end)-1) + Pes_ch_V2G(Ind, 1:Tadd(end)-1)*eff_es - Pes_dis_V2G(Ind, 1:Tadd(end)-1)/eff_es];
                end
                Constraints = [Constraints; Ees_V2G(Ind, Tadd(end)) == Eout_V2G(Ind)];
            end
        end

        % ---------- Voltage Drop and Equipment Limits ----------
        Constraints = [Constraints; Aim * Usqa == 2 * Pl .* repmat(Dis.Rexi, 1, T) + 2 * Ql .* repmat(Dis.Xexi, 1, T)];

        % Power and voltage operating limits
        for t = 1:T
            for i = 1:Nl
                Constraints = [Constraints; cone([Dis.Ecap(i); Pl(i,t); Ql(i,t)])];
            end
        end
        for t = 1:T
            for i = 1:Ns
                Constraints = [Constraints; cone([Scap; Ps(i,t); Qs(i,t)])];
            end
        end

        Constraints = [Constraints; Vi_min^2 <= Usqa <= Vi_max^2];

        % P2P power transaction constraints
        Constraints = [Constraints; Pp2p_buy >= 0];
        Constraints = [Constraints; Pp2p_sell >= 0];
        Constraints = [Constraints; sum(PnetLD(p2pSet, :) - Pfeed(p2pSet, :), 1) == Pp2p_buy - Pp2p_sell];

        % ---------- SCR (Spare Capacity Ratio) and Asset Value Constraints ----------

        % Line-based SCR and PV (Present Value)
        Constraints = [Constraints; -repmat(VirCapL, 1, T) <= Pl <= repmat(VirCapL, 1, T)];
        Constraints = [Constraints; VirCapL >= 0];
        Constraints = [Constraints; SCRL == (Dis.Ecap - VirCapL) ./ Dis.Ecap];
        Constraints = [Constraints; SCRL == sum(delMarL, 2)];
        Constraints = [Constraints; PVL == AssetL + sum(slopeL .* delMarL, 2)];
        Constraints = [Constraints; 0 <= delMarL <= interval];

        % Substation-based SCR and PV (Present Value)
        Constraints = [Constraints; -repmat(VirCapS, 1, T) <= Ps <= repmat(VirCapS, 1, T)];
        Constraints = [Constraints; VirCapS >= 0];
        Constraints = [Constraints; SCRS == (Scap - VirCapS) ./ Scap];
        Constraints = [Constraints; SCRS == sum(delMarS, 2)];
        Constraints = [Constraints; PVS == AssetS + sum(slopeS .* delMarS, 2)];
        Constraints = [Constraints; 0 <= delMarS <= interval];


        %% ---------- Second-Stage Objective and Optimization ----------
        % Redefine objective components for stage 2
        Fdr = sum(sum(delpi0)) * DR_cost;
        Fes = sum(sum(k_degra / 100 * (Pes_dis + Pes_ch) .* ESS_capcost_perunit .* repmat(EShour_duration, 1, T)));
        Fes_V2G = sum(sum(k_degra / 100 * (Pes_dis_V2G + Pes_ch_V2G) .* ESS_capcost_perunit_V2G .* repmat(EShour_duration_V2G, 1, T)));
        Fpdn1 = sum(rou_ToU .* sum(PnetLD(singleSet, :), 1)) - sum(rou_FiT .* sum(Pfeed(singleSet, :), 1));
        Fpdn2 = sum(rou_ToU .* Pp2p_buy) - sum(rou_FiT .* Pp2p_sell);

        % Constraint to ensure second-stage cost does not exceed first-stage by more than permitted deviation
        Constraints = [Constraints; Fes + Fes_V2G + Fpdn1 + Fpdn2 + Fdr <= Fpdn_init * (1 + Permit_costDeviation)];

        % If no P2P trading, fix control variables to first-stage values for singleSet buses
        if p2pRate <= 0
            Constraints = [Constraints; pi0Actu(singleSet, :) == pi0Actu_init(singleSet, :)];
            Constraints = [Constraints; Pes(singleSet, :) == Pes_init(singleSet, :)];
            Constraints = [Constraints; Pes_V2G(singleSet, :) == Pes_V2G_init(singleSet, :)];
            Constraints = [Constraints; PV_out(singleSet, :) == PV_out_init(singleSet, :)];
        end

        % New objective: minimize the present value of infrastructure (line + substation)
        PVofAsset = sum(PVL) + PVS;
        Objective = PVofAsset;

        % Solve the second-stage optimization
        options = sdpsettings('solve', 'mosek');
        options = sdpsettings(options, 'verbose', 2);
        solvesdp(Constraints, Objective, options);

        %% ---------- Output and Save Results ----------

        % Convert costs and variables to numeric values
        Fpdn = double(Fdr + Fes + Fes_V2G + Fpdn1 + Fpdn2);
        Fpdn1 = double(Fpdn1);
        Fpdn2 = double(Fpdn2);
        Fes = double(Fes);
        Fes_V2G = double(Fes_V2G);
        Fdr = double(Fdr);
        Fpdn2Hourly = double(rou_ToU .* Pp2p_buy - rou_FiT .* Pp2p_sell);

        % Power flow and voltages
        Ps = double(Ps);
        Pl = double(Pl);
        Ql = double(Ql);
        Sl = sqrt(Pl.^2 + Ql.^2);

        Usqa = double(Usqa);
        Qs = double(Qs);
        PV_out = double(PV_out);
        PV_cur = PV_pre - PV_out;

        % Energy storage (ESS and V2G)
        Ees = double(Ees);
        Pes = double(Pes);
        Pes_dis = double(Pes_dis);
        Pes_ch = double(Pes_ch);

        Pes_ch_V2G = double(Pes_ch_V2G);
        Pes_dis_V2G = double(Pes_dis_V2G);
        Pes_V2G = double(Pes_V2G);
        Ees_V2G = double(Ees_V2G);

        % P2P market results
        Pp2p_buy = double(Pp2p_buy);
        Pp2p_sell = double(Pp2p_sell);
        PnetLD = double(PnetLD);
        Pfeed = double(Pfeed);

        % Demand response results
        pi0Actu = double(pi0Actu);
        DRrate = sum(pi0Actu, 1) ./ sum(pi0, 1);

        % Final asset value results
        VirCapL = double(VirCapL);
        SCRL = double(SCRL);
        PVL = double(PVL);
        VirCapS = double(VirCapS);
        SCRS = double(SCRS);
        PVS = double(PVS);
        PVA_init;
        PVA = double(PVofAsset);

        toc;
    end
end


