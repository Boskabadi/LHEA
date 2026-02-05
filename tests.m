% for A class sensors use simulink model sensor A, and change XXX in file and simulink

clc
clear all

%% sensor1 and actuator1 parameter
use2ndFcn_X=1; % select 2nd order (1) or 8th order (0) transfer function

min_X = 0; %lower measurement limit
max_X = 5100; %upper measurement limit

%2nd order
T90_X_2 = 5; %response time minutes
T_X_2 = T90_X_2/(24*60)/3.89;

%8th order
T90_X_8 = 20; %response time minutes
T_X_8 = T90_X_8/(24*60)/11.7724;

std_X = 0.025; %standard deviation of noise
noiseseed_X = use2ndFcn_X; %noise seed for random generator (mean=0, std=1, sample=1 per minute)
usenoise_X = 1; %select noise or not (0=no noise, 1=use noise) for sensor
useideal_X = 0; %select ideal sensor or not (0=non-ideal, 1=ideal) for sensor
useidealac_X = 1; %select ideal actuator or not (0=non-ideal, 1=ideal) for sensor

