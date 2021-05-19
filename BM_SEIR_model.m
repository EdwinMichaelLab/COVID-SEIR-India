%% BM SEIR model

%% STEP 1: Define uniform parameter prior min/max values
parameter_priors = ... 
    [0.125 2.0 % 1: beta1
    2.0 2.0 % 2: alpha - set to quickly move S --> R1 at lockdown start
    0.16 0.5 % 3: sigma inverse of 2-6 days
    0.25 0.5 % 4: rho 25-50% of cases are asymptomatic
    0.125 0.33 % 5: gammaA inverse of 4-14 days recovery
    0.125 0.33 % 6: gammaM
    0.125 0.33 % 7: gammaH
    0.125 0.33 % 8: gammaC
    0.05 0.20 % 9: delta1 inverse of 1-10 days - modified to 1-20 days
    0.06 0.25 % 10: delta2 inverse of 4-15 days
    0.09 1 % 11: delta3 inverse of 1-11 days
    0.08 0.25 % 12: m
    3 5 % 13: REPLACED for potential rebase %% lockdown ratio, alpha/lambda 
    0.1 0.3 % 14: epsilon, proportion of symptomatic cases undetected
    0.05 0.3 % 15: x1
    0.2 0.3 % 16: x2 
    0.2 0.8  % 17: x3
    0.25 0.9 % 18: d 1-[0.58 0.85]
    2 500 % 19: E0
    0.5 1 % 20: z
    0.125 0.33 % 21: gammaQ
    0.06 0.25 % 22: deltaQ
    0.05 0.3 % 23: p
    0.1428 2.5]; % 24: beta2

%other parameters
%fit to section of fit from start?
fit_to_section = false;

%number of samples to retain after RMSD sorting
%NOTE: nSelectFitting >= nSelectSimulating
nSelectFitting = 500;
nSelectSimulating = 500;

%make sure >
nSelectFitting = max(nSelectFitting, nSelectSimulating);

%% STEP 2: Randomly sample parameter sets from prior distributions
ParamSets = SampleParamSets(nDraws,parameter_priors);
%save originals for blending
ParamSets_orig = ParamSets;

%% STEP 3: Run the model using sampled parameter sets up to the last time for which we have data
% setup IC
% Set initial compartment values
IA0 = zeros(1,nDraws) + 0/NPop;
IP0 = Confirmed(1) ./ ((1-ParamSets(14,:)) .* ParamSets(9, :) .* NPop);
IM0 = zeros(1,nDraws) + 0/NPop;
IH0 = zeros(1,nDraws) + 0/NPop;
IC0 = zeros(1,nDraws) + 0/NPop;
D0  = zeros(1,nDraws) + Deaths(1)/NPop;
R10 = zeros(1,nDraws) + 0/NPop;
R20 = zeros(1,nDraws) + 0/NPop;
V0 = zeros(1,nDraws) + 0/NPop;
R30 = zeros(1,nDraws) + 0/NPop;
B0 = zeros(1,nDraws) + 0/NPop;

% 2nd Variant
E20 = zeros(1,nDraws) + 0/NPop;
IA20 = zeros(1,nDraws) + 0/NPop;
IP20 = zeros(1,nDraws) + 0/NPop;
IM20 = zeros(1,nDraws) + 0/NPop;
IH20 = zeros(1,nDraws) + 0/NPop;
IC20 = zeros(1,nDraws) + 0/NPop;
R40 = zeros(1,nDraws) + 0/NPop;
D20 = zeros(1,nDraws) + 0/NPop;


% % E0 and S0 defined in parfor loop
INPop = 1/NPop;
E0 = ParamSets(19,:) .* INPop; % E0 = P(18)/NPop;
E20 = ParamSets(19,:) .* INPop;
S0 = (1 - IA0(1) - IP0(1) - IM0(1) - IH0(1) - IC0(1) - D0(1) - R10(1) - R20(1) - R30(1) - V0(1) - B0(1)) - E0 - E20; % S0 = 1 - sum([E0,IA0,IP0,IM0,IH0,IC0,D0,R10,R20]);

% Initalize vaccine refusal class, S2
S20 = S0.*vac_refusal;
S0 = S0 - S20;
%initialize traj. storage
S = S0; S2 = S20; E = E0; IA = IA0; IP = IP0; IM = IM0; IH = IH0;
IC = IC0; D = D0; R1 = R10; R2 = R20; V = V0;
R3 = R30; B = B0;

E2 = E20; IA2 = IA20; IP2 = IP20; IM2 = IM20; IH2 = IH20;
IC2 = IC20;R4 = R40;D2 = D20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup segments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup segments
%equi-spaced?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
segment_length = 10;
loops = floor(length(timeRef) / segment_length);

%lookup for segments
segment_steps = zeros(1,loops + 1);

current_step = 1;

for i = 1:loops
    segment_steps(i) = current_step;
    current_step = current_step + segment_length;
end

%get ending point
if length(timeRef) - current_step < segment_length
    current_step = length(timeRef);
end

segment_steps(end) = current_step;

%user designated?%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%segment_steps = [1 14 length(timeRef)];

%show steps
segment_steps

% get number of segments
num_segments = length(segment_steps);

% and progress steps per segment
prog_steps = floor(50/num_segments);

%loop over segments
for i = 1:(length(segment_steps) - 1)
    current_step = segment_steps(i)
    segment_end = segment_steps(i+1);
    
    % solve first fortnight
    [S0,E0,IA0,IP0,IM0,IH0,IC0,D0,R10,R20,V0,S20,R30,B0,E20,IA20,IP20,IM20,IH20,IC20,R40,D20] = SEIR_covid_model(ParamSets,NPop,...
        S(end,:),E(end,:),IA(end,:),IP(end,:),IM(end,:),IH(end,:),IC(end,:),D(end,:),R1(end,:),R2(end,:),...
        V(end,:),S2(end,:),R3(end,:),B(end,:),E2(end,:),IA2(end,:),IP2(end,:),IM2(end,:), ...
	IH2(end,:),IC2(end,:),R4(end,:),D2(end,:),current_step,segment_end,q,quarantine_start,vac_start,movement_range,prog_flag,M,mr,vacdata,prog_steps);
    %size(S0)

    % prepend IC, the first point is at t0
    S = [S; S0(2:end,:)];
    E = [E; E0(2:end,:)];
    IA = [IA; IA0(2:end,:)];
    IP = [IP; IP0(2:end,:)];
    IM = [IM; IM0(2:end,:)];
    IH = [IH; IH0(2:end,:)];
    IC = [IC; IC0(2:end,:)];
    D = [D; D0(2:end,:)];
    R1 = [R1; R10(2:end,:)];
    R2 = [R2; R20(2:end,:)];
    V = [V; V0(2:end,:)];
    S2 = [S2; S20(2:end,:)];
    R3 = [R3; R30(2:end,:)];
    B = [B; B0(2:end,:)];
    E2 = [E2; E20(2:end,:)];
    IA2 = [IA2; IA20(2:end,:)];
    IP2 = [IP2; IP20(2:end,:)];
    IM2 = [IM2; IM20(2:end,:)];
    IH2 = [IH2; IH20(2:end,:)];
    IC2 = [IC2; IC20(2:end,:)];
    R4 = [R4; R40(2:end,:)];
    D2 = [D2; D20(2:end,:)];
    % Correct S, V for vaccination.
    totalv = 0;
    totalb = 0;
    for vac_day = current_step:segment_end
        if vac_day <= length(Vaccinated)
            totalv = totalv + Vaccinated(vac_day);
        else
            totalv = totalv + mean(Vaccinated(end-21:end));
        end
    end
    
    % S --> V transition
    S(segment_end, :) = S(segment_end, :) - totalv;
    V(segment_end, :) = V(segment_end, :) + totalv;

    % lets grab the reasonable trajectories,
    % we use a form of relative RMSE as a distance metric to indicate how far 
    % the model outputs are from the observed case data. The RMSE is calculated 
    % based on cumulative confirmed case counts and deaths.
    % First, we need to calculate the model-predicted cumulative confirmed 
    % cases. We assume the number of new cases that will be detected each day to 
    % be a fraction of those leaving presymptomatic class and entering the 
    % symptomatic pipeline. We assume epsilon of those are mildly symptomatic 
    % cases that will not get tested/reported.
    epsilon = ParamSets(14,:);
    delta1 = ParamSets(9, :);
    x3 = ParamSets(17, :);
    m = ParamSets(12, :);

    %fit to this section? or from start?
    %if from start we need the pred_C from start
    if fit_to_section
        %pred_C = cumsum((1-epsilon).*IP0.*delta1*NPop);
        pred_C = (1-epsilon).*IP0.*delta1*NPop;
        pred_D = x3.*m.*(IC+IC2)*NPop;
        fit_start = current_step;
    else
        %pred_C = cumsum((1-epsilon).*IP.*delta1*NPop);
        pred_C = (1-epsilon).*IP.*delta1*NPop;
        pred_D = x3.*m.*(IC+IC2)*NPop;
        fit_start = 1;
    end
    
    % B199 Cases RMSE
    pred_C2 = (1-epsilon).*IP2.*delta1*NPop;

    %NOTE: I was fitting to the current segment, seems to work better from
    %the start to the current point!
    % Reformat model predictions and observed data for confirmed
    % cases and deaths
    

    % B199 Cases x1, y1
     x1 = pred_C(1:length(diff(Confirmed(fit_start:segment_end))),:); % model-predicted cumulative cases
     y1 = diff(Confirmed(fit_start:segment_end))'; % observed cumulative cases
     x2 = pred_D(1:length(diff(Confirmed(fit_start:segment_end))),:); % model-predicted deaths
     y2 = diff(Deaths(fit_start:segment_end))'; % observed deaths
     x3 = pred_C2(1:length(diff(Confirmed(fit_start:segment_end))),:); % model-predicted cumulative cases
     y3 = diff(Confirmed(fit_start:segment_end))'; % observed cumulative cases
     RMSE = sqrt(mean([((x1-y1).^2)./std(x1);((x2-y2).^2)./std(x2)]));

    % calculate combined RMSE metric, MSE normalized by standard deviation to
    % avoid higher weighting of confirmed cases than deaths (there are much 
    % fewer deaths than cases!)
    %RMSE = sqrt(mean(((x1-y1).^2)./std(x1))); % Only fit to Confirmed cases for now
    [x,ind] = sort(RMSE);
    
    % select the best nSelectFitting (200) parameter sets based on the corresponding RMSE metric
    id = ind(1:nSelectFitting);
    
    %ParamSets = ParamSets(:,id);
    %Now take selected trajectories 
    % Set initial compartment values
    S = S(:,id); E = E(:,id); IA = IA(:,id); IP = IP(:,id); IM = IM(:,id); 
    IH = IH(:,id); IC = IC(:,id); D = D(:,id); R1 = R1(:,id); R2 = R2(:,id);
    V = V(:,id); S2 = S2(:,id); R3 = R3(:,id);
    B = B(:, id); E2 = E2(:,id); IA2 = IA2(:,id); IP2 = IP2(:,id); IM2 = IM2(:,id); 
    IH2 = IH2(:, id); IC2 = IC2(:,id); R4 = R4(:,id);D2 = D2(:,id);
    
    %and parameters
    ParamSets = ParamSets(:,id);
    
    %if last loop don't bother with replicating traj. and new priors
    if i < length(segment_steps) - 1

        %replicate traj. to get nDraw simulations.
        %number of replications
        nReplicas = nDraws / nSelectFitting;

        S = repmat(S,1,nReplicas);
        E = repmat(E,1,nReplicas);
        IA = repmat(IA,1,nReplicas);
        IP = repmat(IP,1,nReplicas);
        IM = repmat(IM,1,nReplicas);
        IH = repmat(IH,1,nReplicas);
        IC = repmat(IC,1,nReplicas);
        D = repmat(D,1,nReplicas);
        R1 = repmat(R1,1,nReplicas);
        R2 = repmat(R2,1,nReplicas);
        V = repmat(V,1,nReplicas);
        S2 = repmat(S2,1,nReplicas);
        R3 = repmat(R3,1,nReplicas);
        B = repmat(B,1,nReplicas);
        E2 = repmat(E2,1,nReplicas);
        IA2 = repmat(IA2,1,nReplicas);
        IP2 = repmat(IP2,1,nReplicas);
        IM2 = repmat(IM2,1,nReplicas);
        IH2 = repmat(IH2,1,nReplicas);
        IC2 = repmat(IC2,1,nReplicas);
        R4 = repmat(R4,1,nReplicas);
        D2 = repmat(D2,1,nReplicas);


     	% Save parameters?
     	save(sprintf("%i%s.mat", i, Location_arr(1)), "ParamSets");
        %Create new priors? Only if not end of loop
        disp("creating new parameters")
        %get range from selected priors
        %note we have culled ParamSets by here
        %selectedParamSets = ParamSets(:,id);
        parameter_priors_new=[min(ParamSets,[],2) max(ParamSets,[],2)];
        %Randomly sample parameter sets from prior distributions
        ParamSets = SampleParamSets(nDraws,parameter_priors_new);

        %blend with original set
        blend_probability = 0.50;

        %throw rands, and take that proportion from original
        %currently does NOT blend within parameter sets!
        R = rand(1,nDraws);
        ParamSets(:,R<blend_probability) = ParamSets_orig(:,R<blend_probability);
    end
end

disp("End of fitting")

%% STEP 4: Select best-fitting models based on distance metric 
%NOTE: we only need to do this if nSelectFitting > nSelectSimulating

% we use a form of relative RMSE as a distance metric to indicate how far 
% the model outputs are from the observed case data. The RMSE is calculated 
% based on cumulative confirmed case counts and deaths.

% First, we need to calculate the model-predicted cumulative confirmed 
% cases. We assume the number of new cases that will be detected each day to 
% be a fraction of those leaving presymptomatic class and entering the 
% symptomatic pipeline. We assume epsilon of those are mildly symptomatic 
% cases that will not get tested/reported.
epsilon = ParamSets(14,:);
delta1 = ParamSets(9, :);
x3 = ParamSets(17, :);
m = ParamSets(12, :);

pred_C = (1-epsilon).*IP.*delta1*NPop;
pred_D = x3.*m.*IC*NPop;
pred_C2 = (1-epsilon).*IP2.*delta1*NPop;

%pred_C = cumsum((1-epsilon).*IP.*delta1*NPop);

% Define the last time point of data to use for fitting (3 options)

% Option 1: use all data (normal case)
%endt = length(Confirmed); % all data
endt = length(diff(Confirmed)); % all data

x1 = pred_C(1:length(diff(Confirmed(1:endt))),:); % model-predicted cumulative cases
y1 = diff(Confirmed(1:endt))'; % observed cumulative cases
x2 = pred_D(1:length(diff(Confirmed(1:endt))),:); % model-predicted deaths
y2 = diff(Deaths(1:endt))'; % observed deaths
x3 = pred_C2(1:length(diff(Confirmed(1:endt))),:); % model-predicted cumulative cases
y3 = diff(Confirmed(1:endt))'; % observed cumulative cases
     
% calculate combined RMSE metric, MSE normalized by standard deviation to
% avoid higher weighting of confirmed cases than deaths (there are much 
% fewer deaths than cases!)
RMSE = sqrt(mean([((x1-y1).^2)./std(x1);((x2-y2).^2)./std(x2);((x3-y3).^2)./std(x3);]));
%RMSE = sqrt(mean(((x1-y1).^2)./std(x1))); % Only fit to Confirmed cases for now

% select the best nSelectSimulating (200) parameter sets based on the corresponding RMSE metric
[x,i] = sort(RMSE);
id = i(1:nSelectSimulating);
ParamSets = ParamSets(:,id);

% Display median d over the last segment + lockdown ratio on last day
% latest_d = median(ParamSets(18, :)); 
% latest_lockdown = movement_data(length(timeRef));
% if lockdown_flag
%     mr = @(t) lockdown_mod/(1-lockdown_mod);
% end
% 
% if soc_dist_flag
%     ParamSets(18, :) = soc_dist_mod;
% end

%% STEP 5: Run the model for the full simulation period using the sampled parameter sets
% setup IC
% Set initial compartment values
S0 = S(:,id);
E0 = E(:,id);
IA0 = IA(:,id);
IP0 = IP(:,id);
IM0 = IM(:,id);
IH0 = IH(:,id);
IC0 = IC(:,id);
D0 = D(:,id);
R10 = R1(:,id);
R20 = R2(:,id);
V0 = V(:,id);
S20 = S2(:,id);
R30 = R3(:,id);
B0 = B(:,id);
E20 = E2(:,id);
IA20 = IA2(:,id);
IP20 = IP2(:,id);
IM20 = IM2(:,id);
IH20 = IH2(:,id);
IC20 = IC2(:,id);
R40 = R4(:,id);
D20 = D2(:,id);

save(sprintf("final_%s.mat", Location_arr(1)));
%load("final.mat");
% mr = @(t) 0.27/(1-0.27);
% 

for day_segment=length(timeRef):MaxTime-2
    S0 = S(:,id);
    E0 = E(:,id);
    IA0 = IA(:,id);
    IP0 = IP(:,id);
    IM0 = IM(:,id);
    IH0 = IH(:,id);
    IC0 = IC(:,id);
    D0 = D(:,id);
    R10 = R1(:,id);
    R20 = R2(:,id);
    V0 = V(:,id);
    S20 = S2(:,id);
    R30 = R3(:,id);
    B0 = B(:,id);
    E20 = E2(:,id);
    IA20 = IA2(:,id);
    IP20 = IP2(:,id);
    IM20 = IM2(:,id);
    IH20 = IH2(:,id);
    IC20 = IC2(:,id);
    R40 = R4(:,id);
    D20 = D2(:,id);

        
    %solve, start from t0 so t dependent values aligned
       % solve first fortnight
    [S,E,IA,IP,IM,IH,IC,D,R1,R2,V,S2,R3,B,E2,IA2,IP2,IM2,IH2,IC2,R4,D2] = SEIR_covid_model(ParamSets,NPop,...
        S0(end,:),E0(end,:),IA0(end,:),IP0(end,:),IM0(end,:),IH0(end,:),IC0(end,:),D0(end,:),R10(end,:),R20(end,:),...
        V0(end,:),S20(end,:),R30(end,:),B0(end,:),E20(end,:),IA20(end,:),IP20(end,:),IM20(end,:), ...
	IH20(end,:),IC20(end,:),R40(end,:),D20(end,:),day_segment,day_segment+2,q,quarantine_start,vac_start,movement_range,prog_flag,M,mr,vacdata,prog_steps);
    
       % prepend IC, the first point is at t0
    S = [S0; S(2:end-1,:)];
    E = [E0; E(2:end-1,:)];
    IA = [IA0; IA(2:end-1,:)];
    IP = [IP0; IP(2:end-1,:)];
    IM = [IM0; IM(2:end-1,:)];
    IH = [IH0; IH(2:end-1,:)];
    IC = [IC0; IC(2:end-1,:)];
    D = [D0; D(2:end-1,:)];
    R1 = [R10; R1(2:end-1,:)];
    R2 = [R20; R2(2:end-1,:)];
    V = [V0; V(2:end-1,:)];
    S2 = [S20; S2(2:end-1,:)];
    R3 = [R30; R3(2:end-1,:)];
    B = [B0; B(2:end-1,:)];
    E2 = [E20; E2(2:end-1,:)];
    IA2 = [IA20; IA2(2:end-1,:)];
    IP2 = [IP20; IP2(2:end-1,:)];
    IM2 = [IM20; IM2(2:end-1,:)];
    IH2 = [IH20; IH2(2:end-1,:)];
    IC2 = [IC20; IC2(2:end-1,:)];
    R4 = [R40; R4(2:end-1,:)];
    D2 = [D20; D2(2:end-1,:)];  
    
    % Correct S, V for vaccination
    totalv = mean(Vaccinated(end-21:end));

    for paramset=1:500    
        if S(end, paramset) - totalv >= 0
            S(end, paramset) = S(end, paramset) - totalv;
            V(end, paramset) = V(end, paramset) + totalv;
        else
            V(end, paramset) = V(end, paramset) + S(end, paramset);
            S(end, paramset) = 0;
        end

    end
end

% recalcuate the predicted cumulative confirmed cases for full simulation
% period
epsilon = ParamSets(14,:);
delta1 = ParamSets(9, :);
x3 = ParamSets(17, :);
m = ParamSets(12, :);

pred_D = x3.*m.*(IC+IC2)*NPop;
pred_C = (1-epsilon).*IP.*delta1*NPop;
pred_C2 = (1-epsilon).*IP2.*delta1*NPop;
