%% Main.m:

% Updated code can be downloaded from https://github.com/EdwinMichaelLab/covid_model_SEIR.git

% FUNCTION NAME:
%   Main
%
% DESCRIPTION:
%   This is the main function for the project. Calls necessary functions to accept inputs from user, run model, and produce output.

clearvars;
close all
clc;

%% Setup random
rng('default');
rng(1234);

tic
%% INPUTS AND SETUP
%maximum cores to use
M = 128;

% Define counties to model and read data from Johns Hopkins
% Format for each county: "County", "State", "US";
Location_arr = ["India"];
%NPop = 328000000; % 328 million people
NPop = 1366000000;
data = csvread("India.csv");

% Remove today's data because it's a running total and probably low.
data(1, :) = [];

% Find index of last element
last_element = find(data(:, 1) == 0)
last_element = last_element(1);
Confirmed = flip(data(1:last_element, 1));
Deaths = flip(data(1:last_element, 2));
Vaccinated = flip(data(1:last_element, 3));
Vaccinated = [0; diff(Vaccinated)/NPop];
todays_date = datetime('today');
timeRef = (todays_date-last_element):todays_date-1;
Confirmed = Confirmed';
Deaths = Deaths';
Vaccinated = Vaccinated';

% Choose the number of days to simulate beyond latest data (data typically 
% available through yesterday's date)
sim_time = days(datetime(2021, 12, 31) - datetime("Today"));

% Start modeling local epidemic when at least 10 cases were confirmed, effectively
% remove leading zeros from dataset
minNum = min(Confirmed(length(timeRef)-2),10);
Deaths(Confirmed<minNum)=[];
timeRef(Confirmed<minNum)= [];
Vaccinated(Confirmed<minNum)=[];
Confirmed(Confirmed<minNum)=[];
MaxTime = length(timeRef)+sim_time;

days_to_chop = 0;
oldConfirmed = Confirmed;
oldDeaths = Deaths;
save("oldConfirmed.mat")
Deaths = Deaths(1:(end-days_to_chop));
timeRef = timeRef(1:(end-days_to_chop));
Confirmed = Confirmed(1:(end-days_to_chop));
Vaccinated = Vaccinated(1:(end-days_to_chop));
MaxTime = length(timeRef) + sim_time;


%% SMOOOOOOOOOOOOOOOTHER
smoothed_cases = movmean(diff(Confirmed), 7);
smoothed_deaths = movmean(diff(Deaths), 7);

Confirmed = [0 cumsum(smoothed_cases)];
Deaths = [0 cumsum(smoothed_deaths)];

% Lockdown: if flag == 0, maintain lockdown ratio as estimated by today's movement data.
% if 1, change lockdown ratio to lockdown_mod.
lockdown_flag = 0;
lockdown_mod = 0.35;

% Social Distancing: if flag == 0, use "d" parameter as estimated from fit.
% if 1, change d to soc_dist_mod%.
soc_dist_flag = 0; 
soc_dist_mod = 0.30;

% Set quarantine option
% q represents the proportion of asymptomatic/presymptomatic/mild cases that 
% are detected through contact tracing and quarantined, effectively
% preventing them from exposing others and producing new cases
% scenarios: 0, 0.25, 0.5, 0.75
quarantine_start_date = datetime(2020, 08, 27); 
quarantine_start = days(quarantine_start_date - timeRef(1));
q = 0; 

% Vaccination (check vac_rate and vac_efficacy in diff_eqn1.m)
vac_start_date = datetime(2021, 01, 01);
vac_start = days(vac_start_date - timeRef(1));
vac_refusal = 0; % 50% antivax


movement_range = linspace(1, 2, 50000);
mr = [0];
%% Linear interpolation between days
%mr = griddedInterpolant(1:length(movement_ratio_data), movement_ratio_data);

% number of prior parameter sets to sample, default 50000
nDraws = 50000;

% do not allow progress output for web
prog_flag = false;
%load_vac;
%load_2nd_variant;
% ending_timeref = timeRef(end)+sim_time;
% newtimeRef = timeRef(1):ending_timeref;
% [~, vacdates, ~] = intersect(newtimeRef, vactimes);% Produce vacdata
% vacdata = [vacdates vac_v./NPop vac_b./NPop];
vacdata = [0]
%% CALL BM SEIR MODEL
BM_SEIR_model()
toc
%% PLOT OUTPUTS
save(sprintf("%s.mat",Location_arr(1)))
