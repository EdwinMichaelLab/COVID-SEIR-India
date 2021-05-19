%% SEIR_covid_model.m:

% FUNCTION NAME:
%   SEIR_covid_model
%
% DESCRIPTION:
%   This is a helper function for getting the model ready for ode45().
%
% INPUTS:
%   ParamSets: Array of sampled parameters.
%   NPop: Integer, population size
%   MaxTime: Integer, how far to integrate
%   lockdown: Integer, how long to be under lockdown.
%
% OUTPUT:
%   Arrays containing the class values as a function of time.


function [S_out,E_out,IA_out,IP_out,IM_out,IH_out,IC_out,D_out,R1_out,...
		R2_out,V_out,S2_out,R3_out,B_out,E2_out,IA2_out,IP2_out,IM2_out,IH2_out,IC2_out,R4_out,D2_out]...
    = SEIR_covid_model(ParamSets,NPop,...
    S0,E0,IA0,IP0,IM0,IH0,IC0,D0,R10,R20,...
    V0,S20,R30,B0,E20,IA20,IP20,IM20,IH20,IC20,R40,D20,...
    StartTime,MaxTime,q,quarantine_start,vac_start,movement_range,prog_flag,M,mr,vacdata,prog_steps)

%% Initialize simulation

% Set initial compartment values, now passed to integrator

% initialize output arrays
S_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
E_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IA_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IP_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IM_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IH_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IC_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
D_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
R1_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
R2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
V_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
S2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
R3_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
B_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
E2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IA2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IP2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IM2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IH2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
IC2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
R4_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));
D2_out = zeros(MaxTime - StartTime + 1,length(ParamSets(1,:)));

% calculate value to scale progress counts 
pstep = floor(length(ParamSets(1,:)) / prog_steps);

% timing
simu_clk = zeros(1,length(ParamSets(1,:)));

%% Loop through calculations for each parameter set
% can run in parallel for faster computation
parfor (i = 1:length(ParamSets(1,:)), M)
%for i = 1:length(ParamSets(1,:))
    % time loop start
    simu_clk_strt = tic;
    
    %printout progress periodically to 
    %tell web service how far we have progressed
    if prog_flag & mod(i,pstep) == 0
        disp('progress')
    end

    % pull out one set of parameters
    P = ParamSets(:,i); 
    
    
    % numerically integrate the differential equations
    options = odeset('RelTol', 1e-5);
    [t, pop] = ode45(@diff_eqn1,StartTime:1:MaxTime,...
        [S0(i),E0(i),IA0(i),IP0(i),IM0(i),IH0(i),IC0(i),D0(i),...
		R10(i),R20(i),V0(i),S20(i),R30(i),B0(i),E20(i),IA20(i),IP20(i),...
		IM20(i),IH20(i),IC20(i),R40(i),D20(i)],...
        options,...
       [P(1:18)',q,quarantine_start,P(20:24)',vac_start,movement_range(i)],mr,vacdata);
   
    % store the predictions for each compartment for each parameter set
    S_out(:,i) = pop(:,1);
    E_out(:,i) = pop(:,2);
    IA_out(:,i) = pop(:,3);
    IP_out(:,i) = pop(:,4);
    IM_out(:,i) = pop(:,5);
    IH_out(:,i) = pop(:,6);
    IC_out(:,i) = pop(:,7);
    D_out(:,i) = pop(:,8);
    R1_out(:,i) = pop(:,9);
    R2_out(:,i) = pop(:,10);
    V_out(:,i) = pop(:,11);
    S2_out(:,i) = pop(:,12);    
    R3_out(:,i) = pop(:,13);
    B_out(:,i)  = pop(:,14);
    E2_out(:,i) = pop(:,15);
    IA2_out(:,i) = pop(:,16);
    IP2_out(:,i) = pop(:,17);
    IM2_out(:,i) = pop(:,18);
    IH2_out(:,i) = pop(:,19);
    IC2_out(:,i) = pop(:,20);
    R4_out(:,i) = pop(:,21);
    D2_out(:,i) = pop(:,22);

    %time loop end
    simu_clk(i) = toc(simu_clk_strt);
    
end

%print compute time
disp("ODE compute time: " + (mean(simu_clk) * length(ParamSets(1,:))) + " seconds, max " + max(simu_clk) + " seconds.");

end
