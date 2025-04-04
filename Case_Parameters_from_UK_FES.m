 % Load network topology and profile data
load(['SSNG_Keydata_network_' num2str(Network_ID)]);
load('PV_typicalProfiles.mat');
load('GenericLoad_ProfileClass.mat');

%% Voltage / Power Flow base parameters
Vbase = 230;    
Sbase = 1000;   
Ibase = Sbase / Vbase;
Zbase = Vbase^2 / Sbase;
Vi_min = 0.93;
Vi_max = 1.1;
Ce = 100 / 1000;          

% Substation settings
Scap = 500;              
Ns = 1;
AssetS = 40500;  

%% Line parameters initialization
Dis.st = Upstream_Busbar_Code;
Dis.ed = Downstream_Busbar_Code;
Dis.Rexi = zeros(length(Network_Cable_Type), 1);
Dis.Xexi = zeros(length(Network_Cable_Type), 1);
Dis.Ecap = zeros(length(Network_Cable_Type), 1);

% Calculate per-line resistance, reactance, and ampacity
for m = 1:length(Network_Cable_Type)
    for typeInd = 1:length(Set_Cable_Type)
        if Network_Cable_Type(m) == Set_Cable_Type(typeInd)
            linetype = typeInd;
            break
        end
    end
    Dis.Ecap(m, 1) = Set_Cable_RXA(linetype, 5) / Ibase;
    Dis.Rexi(m, 1) = Set_Cable_RXA(linetype, 1) * Network_Cable_Length(m) / Zbase;
    Dis.Xexi(m, 1) = Set_Cable_RXA(linetype, 2) * Network_Cable_Length(m) / Zbase;
end

% Assign investment cost per cable type
Price_perMeter = zeros(length(Network_Cable_Type), 1);
for m = 1:length(Network_Cable_Type)
    switch Network_Cable_Type(m)
        case '4c_0.06'
            Price_perMeter(m,1) = 46.7;
        case '4c_0.1'
            Price_perMeter(m,1) = 47.6;
        case '4c_0.2'
            Price_perMeter(m,1) = 50.4;
        case '4c_0.4'
            Price_perMeter(m,1) = 51.4;
    end
end

AssetL = Price_perMeter .* Network_Cable_Length;  % Line investment cost

Nl = length(Dis.st);
if max(Dis.st) > 250
    busSet = unique([5; Dis.ed]);  
else
    busSet = unique([1; Dis.ed]);
end
Nb = length(busSet);

Dis.shorst = zeros(Nl,1);
Dis.shored = zeros(Nl,1);
if max(Dis.st) > 500
    for l = 1:Nl
        Dis.shorst(l) = find(busSet == Dis.st(l));
        Dis.shored(l) = find(busSet == Dis.ed(l));
    end
else
    Dis.shorst = Dis.st;
    Dis.shored = Dis.ed;
end

Aim = zeros(Nl, Nb);
for m = 1:Nl
    Aim(m, Dis.shorst(m)) = 1;
    Aim(m, Dis.shored(m)) = -1;
end


%% PV Allocation and Capacity Setup
noPV_node_init = find(PV_Potential == 0);
NumInc = round(0.4 * length(noPV_node_init));

if III == 0 || (III >= 2001 && III <= 2099)                            
    temphh = noPV_node_init(randperm(length(noPV_node_init), NumInc));
    temphh_init = temphh;
elseif III == -1 || (III >= 1001 && III <= 1099) || (III >= 3001 && III <= 3099)
    temphh = temphh_init;
end
PV_Potential(temphh) = 1;

shorbus_PV = find([0; PV_Potential] > 0.1);
PVnum_eachnode = [0; PV_Potential];

if PVpene_rateSet(PVpene_rateInd) >= 1.0
    PVpene_rate = 1;
else
    PVpene_rate = PVpene_rateSet(PVpene_rateInd);
end

if III == 0 || (III >= 2001 && III <= 2099)                         
    shorbus_PVlocAtoE = shorbus_PV(sort(randperm(length(shorbus_PV), ceil(length(shorbus_PV) * PVpene_rate))));
    shorbus_PVloc_init = shorbus_PVlocAtoE;
elseif III == -1 || (III >= 1001 && III <= 1099) || (III >= 3001 && III <= 3099)
    shorbus_PVlocAtoE = shorbus_PVloc_init;
end

PVnum_eachnode(setdiff([1:Nb]', shorbus_PVlocAtoE)) = 0;

Connected_PV_Aug = [0; Connected_PV];
PVnum_eachnode(find(Connected_PV_Aug >= 0.01)) = Connected_PV_Aug(find(Connected_PV_Aug >= 0.01));

shorbus_PVlocFinal = find(PVnum_eachnode >= 0.1);

if PVpene_rateSet(PVpene_rateInd) >= 1.0
    PV_cap_base = 4.5 * PVpene_rateSet(PVpene_rateInd);
else
    PV_cap_base = 4.2;
end

if III == 0 || (III >= 2001 && III <= 2099)                            
    PV_cap = ones(Nb, 1) * 0.93 * PV_cap_base + 0.14 * PV_cap_base * rand(Nb, 1);
    PV_cap_init = PV_cap;
elseif III == -1 || (III >= 1001 && III <= 1099) || (III >= 3001 && III <= 3099)
    PV_cap = PV_cap_init;
end

if LDsceInd >= 1 && LDsceInd <= 3
    PV_pre = PVnum_eachnode .* PV_cap * PV_Autumn_typical(PVsceInd, :);
elseif LDsceInd >= 4 && LDsceInd <= 6
    PV_pre = PVnum_eachnode .* PV_cap * PV_Highsummer_typical(PVsceInd, :);
elseif LDsceInd >= 7 && LDsceInd <= 9
    PV_pre = PVnum_eachnode .* PV_cap * PV_summer_typical(PVsceInd, :);
elseif LDsceInd >= 10 && LDsceInd <= 12
    PV_pre = PVnum_eachnode .* PV_cap * PV_spring_typical(PVsceInd, :);
elseif LDsceInd >= 13 && LDsceInd <= 15
    PV_pre = PVnum_eachnode .* PV_cap * PV_Winter_typical(PVsceInd, :);
end


%% Load  and profile
% Load scaling based on scenario index (rateIIInd)
if rateIIInd == 1 || rateIIInd == 8
    Cload = [0; Connected_loads] * 0.5296;
elseif rateIIInd == 2
    Cload = [0; Connected_loads] * 0.5182;
elseif rateIIInd == 3
    Cload = [0; Connected_loads] * 0.5826;
elseif rateIIInd == 4
    Cload = [0; Connected_loads] * 0.6783;
elseif rateIIInd == 5
    Cload = [0; Connected_loads] * 0.8274;
elseif rateIIInd == 6
    Cload = [0; Connected_loads] * 0.8766;
elseif rateIIInd == 7
    Cload = [0; Connected_loads] * 0.8774;
elseif rateIIInd == 9
    Cload = [0; Connected_loads] * 0.5165;
elseif rateIIInd == 10
    Cload = [0; Connected_loads] * 0.5983;
elseif rateIIInd == 11
    Cload = [0; Connected_loads] * 0.7535;
elseif rateIIInd == 12
    Cload = [0; Connected_loads] * 1.0037;
elseif rateIIInd == 13
    Cload = [0; Connected_loads] * 1.1372;
elseif rateIIInd == 14
    Cload = [0; Connected_loads] * 1.1481;
end

Cload_origin = Cload;
Cload(shorbus_PVlocFinal) = Cload_origin(shorbus_PVlocFinal);  

Rtemp = (sum(Cload_origin) - sum(Cload(shorbus_PVlocFinal))) / (sum(Cload_origin) - sum(Cload_origin(shorbus_PVlocFinal)));
Cload(setdiff(1:Nb, shorbus_PVlocFinal)) = Cload_origin(setdiff(1:Nb, shorbus_PVlocFinal)) * Rtemp;

% load profiles
Max_GenericLoad_ProfileClass1 = max(max(GenericLoad_ProfileClass1));
pu_Profile_GenericLoad_ProfileClass1 = GenericLoad_ProfileClass1 / Max_GenericLoad_ProfileClass1;
pu_Profile_GenericLoad_ProfileClass1(:, 13) = pu_Profile_GenericLoad_ProfileClass1(:, 13);  % Preserve original

pu_Profile_24hour = zeros(24, size(pu_Profile_GenericLoad_ProfileClass1, 2));
for h = 1:24
    pu_Profile_24hour(h, :) = mean(pu_Profile_GenericLoad_ProfileClass1(2*(h-1)+1:2*h, :));
end

% Calculate active and reactive load profiles (pi0, qi0)
pi0 = Cload * pu_Profile_24hour(:, LDsceInd)';
qi0 = pi0 / 0.98 * sqrt(1 - 0.98^2);  % Assume 0.98 power factor


%% Energy Storage System (ESS) Setup
ESS_rate = ESSrateSet(ESSrateInd); 
ESS_num = ceil(length(find(PVnum_eachnode >= 0.1)) * ESS_rate);  
eff_es = 0.97;  

smax = 16 * ones(Nb, 1) + 3.5 * rand(Nb, 1);
cmax = 5 * ones(Nb, 1) + 1.2 * rand(Nb, 1);

% Select ESS installation locations
if ESS_num <= length(shorbus_PVlocFinal)
    if III == 0 || (III >= 2001 && III <= 2099)
        ESS_loc = shorbus_PVlocFinal(sort(randperm(length(shorbus_PVlocFinal), ESS_num)));
        ESS_loc_init = ESS_loc;
    elseif III == -1 || (III >= 1001 && III <= 1099) || (III >= 3001 && III <= 3099)
        ESS_loc = ESS_loc_init;
    end
end

% Zero out capacity where no ESS is installed
smax(setdiff([1:Nb]', ESS_loc)) = 0;
cmax(setdiff([1:Nb]', ESS_loc)) = 0;

% For CT scenarios (ESSrateInd 8–14), capacity is scaled to 87%
if ESSrateInd >= 8 && ESSrateInd <= 14
    cmax = cmax * 0.87;
    smax = smax * 0.87;
end

if (III >= 3001) && (III <= 3099)
    cmax = cmax * 1e-3;
    smax = smax * 1e-3;
end

EShour_duration = zeros(Nb, 1);
EShour_duration(find(cmax >= 0.1)) = smax(find(cmax >= 0.1)) ./ cmax(find(cmax >= 0.1));   

k_degra = 0.01544;  % Degradation factor

ESS_capcost_perunit = 100;  

LDSetweekday = [1, 4, 7, 10, 13];
LDSetweekend = setdiff(1:15, LDSetweekday);


%% Peer-to-Peer (P2P) Trading Set Definition
if III == 0 || (III >= 2001 && III <= 2099) || (III >= 3001 && III <= 3099)
    p2pRate = 0;
elseif III == -1 || (III >= 1001 && III <= 1099)
    p2pRate = 1;
end

% Classify load and PV node sets
loadSet = find(Cload > 0.01);
lowloadSet = find((mean(Cload) > Cload) & (Cload > 0.01));
pvSet = shorbus_PVlocFinal;

singloadSet = setdiff(loadSet, pvSet);
singpvSet = setdiff(pvSet, loadSet);
pv_loadSet = intersect(pvSet, loadSet);

% Nodes for P2P participation
p2pSet_singload = singloadSet(sort(randperm(length(singloadSet), ceil(length(singloadSet) * p2pRate))));
p2pSet_singpv = singpvSet(sort(randperm(length(singpvSet), ceil(length(singpvSet) * p2pRate))));
p2pSet_pv_load = pv_loadSet(sort(randperm(length(pv_loadSet), ceil(length(pv_loadSet) * p2pRate))));
p2pSet = [1; p2pSet_singload; p2pSet_singpv; p2pSet_pv_load];
singleSet = setdiff([1:Nb]', p2pSet);  % Nodes not in P2P

%% Piecewise Linearization for Asset Valuation
discountR = 0.06;
loadR = 0.02;
step = 40;
aax = linspace(0, 1, step+1); interval = aax(2) - aax(1);
aay = AssetL .* (1 + discountR) .^ (log(1 - aax) / log(1 + loadR));
aaz = AssetS * (1 + discountR) .^ (log(1 - aax) / log(1 + loadR));
slopeL = (aay(:,2:end) - aay(:,1:end-1)) / interval;
slopeS = (aaz(:,2:end) - aaz(:,1:end-1)) / interval;

%% Demand Response (DR) Modeling
DR_cost = 0.075;
gap = 0.40;

if (III >= 3001) && (III <= 3099)
    gap = gap * 1e-3;
end

noDRnum = ceil((1 - DR_rateSet(DR_rateInd)) * Nb);
if length(shorbus_PVlocFinal) < noDRnum
    noDRset = shorbus_PVlocFinal;
    noPVset = setdiff([1:Nb]', shorbus_PVlocFinal);
    exnoDR = noPVset(sort(randperm(length(noPVset), noDRnum - length(shorbus_PVlocFinal))));
    noDRset = [noDRset; exnoDR];
else
    noDRset = sort(shorbus_PVlocFinal(randperm(length(shorbus_PVlocFinal), noDRnum)));
end

pi0min = (1 - gap) * pi0; 
pi0max = 2 * repmat(max(pi0, [], 2), 1, T);
pi0min(noDRset, :) = pi0(noDRset, :);  
pi0max(noDRset, :) = pi0(noDRset, :);


%% -----------------------------------EV Charging Load Setup
% Probability distribution of EV charging start times across 24 hours
Proba_Char_stTime = [
    0.0166; 0; 0; 0; 0; 0; 0; 0;
    0.0166; 0.0166; 0.0166; 0.0420; 0.0664; 0.0497; 0.0664; 0.0826;
    0.0993; 0.1112; 0.0918; 0.0750; 0.0664; 0.0828; 0.0500; 0.0500
];

% Number of EVs in each charging category based on penetration scenario
if rateIIInd == 1 || rateIIInd == 8
    EVnum_DumbCharge  = 1;   
    EVnum_SmartCharge = 1;   
    EVnum_V2G         = 1;   
elseif rateIIInd == 2
    EVnum_DumbCharge  = 4;
    EVnum_SmartCharge = 4;
    EVnum_V2G         = 1;
elseif rateIIInd == 3
    EVnum_DumbCharge  = 9;
    EVnum_SmartCharge = 18;
    EVnum_V2G         = 4;
elseif rateIIInd == 4
    EVnum_DumbCharge  = 7;
    EVnum_SmartCharge = 54;
    EVnum_V2G         = 16;
elseif rateIIInd == 5
    EVnum_DumbCharge  = 9;
    EVnum_SmartCharge = 60;
    EVnum_V2G         = 28;
elseif rateIIInd == 6
    EVnum_DumbCharge  = 9;
    EVnum_SmartCharge = 51;
    EVnum_V2G         = 35;
elseif rateIIInd == 7
    EVnum_DumbCharge  = 7;
    EVnum_SmartCharge = 36;
    EVnum_V2G         = 35;
elseif rateIIInd == 9
    EVnum_DumbCharge  = 4;
    EVnum_SmartCharge = 3;
    EVnum_V2G         = 1;
elseif rateIIInd == 10
    EVnum_DumbCharge  = 11;
    EVnum_SmartCharge = 18;
    EVnum_V2G         = 2;
elseif rateIIInd == 11
    EVnum_DumbCharge  = 13;
    EVnum_SmartCharge = 55;
    EVnum_V2G         = 9;
elseif rateIIInd == 12
    EVnum_DumbCharge  = 18;
    EVnum_SmartCharge = 69;
    EVnum_V2G         = 16;
elseif rateIIInd == 13
    EVnum_DumbCharge  = 18;
    EVnum_SmartCharge = 67;
    EVnum_V2G         = 21;
elseif rateIIInd == 14
    EVnum_DumbCharge  = 17;
    EVnum_SmartCharge = 59;
    EVnum_V2G         = 26;
end



%% Dumb Charging of EVs
SimuCoe = 0.30;  
if (III >= 3001) && (III <= 3099)         
    EVnum_DumbCharge = EVnum_DumbCharge + EVnum_SmartCharge + EVnum_V2G;
end

% Dumb EV battery capacity (kWh), normally distributed with clipping
if III == 0 || (III >= 2001 && III <= 2099)                           
    mu = 40; sigma = 2.5;
    dumbEVcap = round(normrnd(mu, sigma, EVnum_DumbCharge, 1));  
    stran = find(dumbEVcap >= mu + 3*sigma | dumbEVcap <= mu - 3*sigma);
    dumbEVcap(stran) = mu;
    dumbEVcap_init = dumbEVcap;
elseif III == -1 || (III >= 1001 && III <= 1099)
    dumbEVcap = dumbEVcap_init;
elseif (III >= 3001) && (III <= 3099)
    mu = 40; sigma = 1.5;
    dumbEVcap = round(normrnd(mu, sigma, EVnum_DumbCharge, 1));  
    stran = find(dumbEVcap >= mu + 3*sigma | dumbEVcap <= mu - 3*sigma);
    dumbEVcap(stran) = mu;
end

% Dumb EV charging power (kW), normally distributed with clipping
if III == 0 || (III >= 2001 && III <= 2099)                            
    mu = 6; sigma = 0.3;
    kW_dumbChar = normrnd(mu, sigma, EVnum_DumbCharge, 1); 
    stran = find(kW_dumbChar >= mu + 3*sigma | kW_dumbChar <= mu - 3*sigma);
    kW_dumbChar(stran) = mu;
    kW_dumbChar_init = kW_dumbChar;
elseif III == -1 || (III >= 1001 && III <= 1099)
    kW_dumbChar = kW_dumbChar_init;
elseif (III >= 3001) && (III <= 3099)
    mu = 6; sigma = 0.3;
    kW_dumbChar = normrnd(mu, sigma, EVnum_DumbCharge, 1); 
    stran = find(kW_dumbChar >= mu + 3*sigma | kW_dumbChar <= mu - 3*sigma);
    kW_dumbChar(stran) = mu;
end

% State of charge (SOC) at start and end of charging
if III == 0 || (III >= 2001 && III <= 2099)                           
    dumbChar_SOC_ST = 0.20 + 0.10 * rand(EVnum_DumbCharge, 1);  
    dumbChar_SOC_ED = 0.93 + 0.07 * rand(EVnum_DumbCharge, 1);  
    dumbChar_SOC_ST_init = dumbChar_SOC_ST;
    dumbChar_SOC_ED_init = dumbChar_SOC_ED;
elseif III == -1 || (III >= 1001 && III <= 1099)
    dumbChar_SOC_ST = dumbChar_SOC_ST_init;
    dumbChar_SOC_ED = dumbChar_SOC_ED_init;
elseif (III >= 3001) && (III <= 3099)
    dumbChar_SOC_ST = 0.20 + 0.10 * rand(EVnum_DumbCharge, 1);  
    dumbChar_SOC_ED = 0.93 + 0.07 * rand(EVnum_DumbCharge, 1);  
end

dumbChar_hours = round((dumbChar_SOC_ED - dumbChar_SOC_ST) .* dumbEVcap ./ kW_dumbChar);

% Charging start time sampled from distribution
if III == 0 || (III >= 2001 && III <= 2099)                           
    dumbChar_stTime = gendist(Proba_Char_stTime', EVnum_DumbCharge, 1);    
    dumbChar_stTime_init = dumbChar_stTime;
elseif III == -1 || (III >= 1001 && III <= 1099)
    dumbChar_stTime = dumbChar_stTime_init;
elseif (III >= 3001) && (III <= 3099)
    dumbChar_stTime = gendist(Proba_Char_stTime', EVnum_DumbCharge, 1); 
end

dumbChar_edTime_Ori = dumbChar_stTime + dumbChar_hours;
evInd_trav = find(dumbChar_edTime_Ori >= 25);  % Vehicles whose charging spans to next day
dumbChar_edTime = dumbChar_edTime_Ori;
dumbChar_edTime(evInd_trav) = 24;

dumbChar_adTime = [];
for i = 1:length(evInd_trav)
    dumbChar_adTime(i).one = 1:1:(dumbChar_edTime_Ori(evInd_trav(i)) - 24);
end

Char_Pload = zeros(EVnum_DumbCharge, 24);
for i = 1:EVnum_DumbCharge
    Char_Pload(i, dumbChar_stTime(i):dumbChar_edTime(i)) = kW_dumbChar(i);
end
for m = 1:length(evInd_trav)
    Char_Pload(evInd_trav(m), dumbChar_adTime(m).one) = kW_dumbChar(m);
end

DumbChar_Pload = Char_Pload * SimuCoe;
if (III >= 3001) && (III <= 3099)
    DumbChar_Pload = Char_Pload * SimuCoe;  
end


%% Smart Charging of EVs
if (III >= 3001) && (III <= 3099)
    EVnum_SmartCharge = 1;
end

% Smart EV battery capacity (kWh)
mu = 40; sigma = 2.5;
if III == 0 || (III >= 2001 && III <= 2099)                           
    smartEVcap = round(normrnd(mu, sigma, EVnum_SmartCharge, 1)); 
    stran = find(smartEVcap >= mu + 3*sigma | smartEVcap <= mu - 3*sigma);
    smartEVcap(stran) = mu;  % Limit capacity within 3σ range
    smartEVcap_init = smartEVcap;
elseif III == -1 || (III >= 1001 && III <= 1099)
    smartEVcap = smartEVcap_init;
elseif (III >= 3001) && (III <= 3099)
    smartEVcap = 5;
end

% Smart EV charging power (kW)
if III == 0 || (III >= 2001 && III <= 2099)                           
    mu = 6; sigma = 0.3;
    kW_smartChar = round(normrnd(mu, sigma, EVnum_SmartCharge, 1)); 
    stran = find(kW_smartChar >= mu + 3*sigma | kW_smartChar <= mu - 3*sigma);
    kW_smartChar(stran) = mu;
    kW_smartChar_init = kW_smartChar;
elseif III == -1 || (III >= 1001 && III <= 1099)
    kW_smartChar = kW_smartChar_init;
elseif (III >= 3001) && (III <= 3099)
    kW_smartChar = 1;
end

% SOC at start and end
if III == 0 || (III >= 2001 && III <= 2099)                           
    smartChar_SOC_ST = 0.20 + 0.10 * rand(EVnum_SmartCharge, 1);
    smartChar_SOC_ED = 0.93 + 0.07 * rand(EVnum_SmartCharge, 1);
    smartChar_SOC_ST_init = smartChar_SOC_ST;
    smartChar_SOC_ED_init = smartChar_SOC_ED;
elseif III == -1 || (III >= 1001 && III <= 1099)
    smartChar_SOC_ST = smartChar_SOC_ST_init;
    smartChar_SOC_ED = smartChar_SOC_ED_init;
elseif (III >= 3001) && (III <= 3099)
    smartChar_SOC_ST = 0.20 + 0.10 * rand(EVnum_SmartCharge, 1);
    smartChar_SOC_ED = 0.93 + 0.07 * rand(EVnum_SmartCharge, 1); 
end

smartChar_hours = round((smartChar_SOC_ED - smartChar_SOC_ST) .* smartEVcap ./ kW_smartChar);

% Charging start time
if III == 0 || (III >= 2001 && III <= 2099)                           
    mu = 1.5; sigma = 5.5;
    smartChar_stTime_ori = round(normrnd(mu, sigma, EVnum_SmartCharge, 1)); 
    stran = find(smartChar_stTime_ori >= mu + 3*sigma | smartChar_stTime_ori <= mu - 3*sigma);
    smartChar_stTime_ori(stran) = mu;
    smartChar_stTime_ori_init = smartChar_stTime_ori;
elseif III == -1 || (III >= 1001 && III <= 1099)
    smartChar_stTime_ori = smartChar_stTime_ori_init;
elseif (III >= 3001) && (III <= 3099)
    mu = 1.5; sigma = 5.5;
    smartChar_stTime_ori = round(normrnd(mu, sigma, EVnum_SmartCharge, 1)); 
    stran = find(smartChar_stTime_ori >= mu + 3*sigma | smartChar_stTime_ori <= mu - 3*sigma);
    smartChar_stTime_ori(stran) = mu;
end

smartChar_stTime = smartChar_stTime_ori;
nega = find(smartChar_stTime_ori <= 0);  
smartChar_stTime(nega) = smartChar_stTime(nega) + 24;

smartChar_edTime_Ori = smartChar_stTime + smartChar_hours;
evInd_trav = find(smartChar_edTime_Ori >= 25);  
smartChar_edTime = smartChar_edTime_Ori;
smartChar_edTime(evInd_trav) = 24;

smartChar_adTime = [];
for i = 1:length(evInd_trav)
    smartChar_adTime(i).one = 1:1:(smartChar_edTime_Ori(evInd_trav(i)) - 24);
end

% charging power profile
Char_Pload = zeros(EVnum_SmartCharge, 24);
for i = 1:EVnum_SmartCharge
    Char_Pload(i, smartChar_stTime(i):smartChar_edTime(i)) = kW_smartChar(i);
end
for m = 1:length(evInd_trav)
    Char_Pload(evInd_trav(m), smartChar_adTime(m).one) = kW_smartChar(evInd_trav(m));
end
SmartChar_Pload = Char_Pload * SimuCoe;



%% V2G Effect of EVs
ESS_capcost_perunit_V2G = 100;

if (III >= 3001) && (III <= 3099)
    EVnum_V2G = 1;
end

% Battery capacity, charge power, SOC and limits
if III == 0 || (III >= 2001 && III <= 2099)                          
    mu = 40; sigma = 2.5;
    smax_V2G = round(normrnd(mu, sigma, EVnum_V2G, 1)); 
    stran = find(smax_V2G >= mu + 3*sigma | smax_V2G <= mu - 3*sigma); 
    smax_V2G(stran) = mu;  
    cmax_V2G = 5 + 2 * rand(EVnum_V2G, 1);
    Ein_V2G = (0.20 + 0.10 * rand(EVnum_V2G, 1)) .* smax_V2G;
    Eout_V2G = (0.93 + 0.07 * rand(EVnum_V2G, 1)) .* smax_V2G;

    smax_V2G_init = smax_V2G;
    cmax_V2G_init = cmax_V2G;
    Ein_V2G_init = Ein_V2G;
    Eout_V2G_init = Eout_V2G;
    
elseif III == -1 || (III >= 1001 && III <= 1099)
    smax_V2G = smax_V2G_init;
    cmax_V2G = cmax_V2G_init;
    Ein_V2G = Ein_V2G_init;
    Eout_V2G = Eout_V2G_init;

elseif (III >= 3001) && (III <= 3099)  
    smax_V2G = 5;
    cmax_V2G = 1;
    Ein_V2G = smax_V2G * 0.25;
    Eout_V2G = smax_V2G * 0.95;
end

EShour_duration_V2G = zeros(EVnum_V2G, 1);
EShour_duration_V2G(find(cmax_V2G >= 0.1)) = smax_V2G(find(cmax_V2G >= 0.1)) ./ cmax_V2G(find(cmax_V2G >= 0.1));

% EV connection time
if III == 0 || (III >= 2001 && III <= 2099)                          
    EV_inTime = gendist(Proba_Char_stTime', EVnum_V2G, 1);
    mu = 10; sigma = 1;
    EV_connecHours = round(normrnd(mu, sigma, EVnum_V2G, 1)); 
    stran = find(EV_connecHours >= mu + 3*sigma | EV_connecHours <= mu - 3*sigma); 
    EV_connecHours(stran) = mu;
    
    EV_inTime_init = EV_inTime;
    EV_connecHours_init = EV_connecHours;

elseif III == -1 || (III >= 1001 && III <= 1099)
    EV_inTime = EV_inTime_init;
    EV_connecHours = EV_connecHours_init;

elseif (III >= 3001) && (III <= 3099)  
    EV_inTime = gendist(Proba_Char_stTime', EVnum_V2G, 1);
    mu = 10; sigma = 1;
    EV_connecHours = round(normrnd(mu, sigma, EVnum_V2G, 1)); 
    stran = find(EV_connecHours >= mu + 3*sigma | EV_connecHours <= mu - 3*sigma); 
    EV_connecHours(stran) = mu;
end

EV_outTime_ori = EV_inTime + EV_connecHours;
evInd_trav = find(EV_outTime_ori >= 25);
EV_outTime = EV_outTime_ori;
EV_outTime(evInd_trav) = 24;

EV_adTime = [];
EV_adTime(EVnum_V2G).one = [];
for i = 1:length(evInd_trav)
    EV_adTime(evInd_trav(i)).one = 1:1:(EV_outTime_ori(evInd_trav(i)) - 24);
end

length_adTime = zeros(EVnum_V2G, 1);
for i = 1:EVnum_V2G
    length_adTime(i) = length(EV_adTime(i).one);
end
evTranver = find(length_adTime > 0);     
evNoTranver = find(length_adTime <= 0);  

% Available time sets for each EV
EV_connecTime = [];
EV_disconnecTime = [];
for i = 1:EVnum_V2G
    EV_connecTime(i).one = union(EV_inTime(i):EV_outTime(i), EV_adTime(i).one);
    EV_disconnecTime(i).one = setdiff([1:24], EV_connecTime(i).one);
end


%% Pev_Dumb_Smart curve: Assign dumb and smart EV charging load to specific buses
Nload = length(lowloadSet);
Pev_Dumb_Smart = zeros(Nb, 24);  

if III == 0 || (III >= 2001 && III <= 2099)                            
    evloadInd = randperm(Nload, EVnum_DumbCharge + EVnum_SmartCharge + EVnum_V2G)';  
    evLoc = lowloadSet(evloadInd);  
    evLoc_init = evLoc;
elseif III == -1 || (III >= 1001 && III <= 1099) || (III >= 3001 && III <= 3099)   
    evLoc = evLoc_init;
end

if (III >= 3001) && (III <= 3099)
    Pev_Dumb_Smart(evLoc(1:end), :) = DumbChar_Pload;
    evLoc_V2G = evLoc(end); 
else 
    Pev_Dumb_Smart(evLoc(1:EVnum_DumbCharge + EVnum_SmartCharge), :) = [DumbChar_Pload; SmartChar_Pload];
    evLoc_V2G = evLoc(EVnum_DumbCharge + EVnum_SmartCharge + 1:end);
end

evLoc_V2G_sort = sort(evLoc_V2G);
Hev = zeros(Nb, EVnum_V2G);
for i = 1:EVnum_V2G
    Hev(evLoc_V2G_sort(i), i) = 1;
end

%% ---------------------------------ToU price
if ToUPri==1
    rou_ToU=zeros(1,T);
    rou_ToU(1:7)=16.44;
    rou_ToU(8:16)=16.44;
    rou_ToU(17:21)=32.55;
    rou_ToU(17)=16.44;
    rou_ToU(22:24)=16.44;
    rou_ToU=rou_ToU/100*1;
elseif ToUPri==0
    rou_ToU=0.2*ones(1,24);
end
rou_ToU=rou_ToU*1.2;
rou_FiT=rou_ToU/4.3166;


