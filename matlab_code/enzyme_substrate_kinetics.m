close all; clear all; clc; 

syms t; 
t_intvl = linspace(0,60,61);
%YP = 1.4*10^(6); % avg Number of cells in colony in YPD
avg_prot_prod = 13000; % avg  number of proteins/cell/sec
cyp_frac = 0.01;
time_span = t*60*60; %seconds for t hours

%% Standard yeast growth
% setting up parameters (parameters values are based on paper Pylvanainen (2005))
% lower-level parameters (by obtaining the physical parameters above from
% experiment, we can solve for the lower-level parameters 
B0 = 4.5; 
B1 = -50;
B2 = 0.3;
B3 = 2; % 2/3: von Bertalanffy function; 2: logistic model 
D = -3; 

B01 = 6;
B02 = 2;

B11 = -30;
B12 = -70;

B21 = 0.9;
B22 = 0.1;

B31 = 6;
B32 = 12;

% champman-Richards model (B0,B2 > 0, B3 > 1, and B1 < 1 - B3 or B0,B2 > 0, 0 < B3 < 1, and 1 - B3 < B1 < 1)
gt = B0*(1 - (B1*exp((-B2)*t)))^(1/(1 - B3)) + D;
gt1 = B0*(1 - (B1*exp((-B21)*t)))^(1/(1 - B31)) + D; %log(Nt);
gt2 = B0*(1 - (B1*exp((-B22)*t)))^(1/(1 - B32)) + D; %log(Nt);
N_t = 10^(gt);
N_t1 = 10^(gt1);
N_t2 = 10^(gt2);
%val_gt = subs(gt, t_intvl);

%% CYP production after a given time_span
CYP = N_t * avg_prot_prod * time_span * cyp_frac; 
CYP1 = N_t1 * avg_prot_prod * time_span * cyp_frac; 
CYP2 = N_t2 * avg_prot_prod * time_span * cyp_frac; 
CYP_t = subs(CYP, t_intvl);




%{
Overview: current version: simple differential rate equations

The implementation is based on https://www.mathworks.com/help/simbio/ug/defining-reaction-rates-with-enzyme-kinetics.html
%}  

%reactions: none
%reaction rate: none
%initial species amount: 
S_i =  8;   %mole
E_i =  CYP2;   %mole
ES_i =  0;  %mole
P_i =  0;   %mole
%starting_amounts = [S_i, E_i, ES_i, P_i];
[t,x] = ode15s(@cal_rates, [0 10], [S_i, subs(E_i,t_intvl), ES_i, P_i]);

%% plot
figure(1);
plot(t,x);
xlabel('time');
ylabel('species amounts');
legend('FC', 'CYP', 'CYP-FC complex', 'FC metabolites');
title('CYP-FC interaction');

%function to calculate the rates of reactions
function [rates] = cal_rates (~, initial_stuff_amounts)
S =  initial_stuff_amounts(1);   %mole
E =  initial_stuff_amounts(2);   %mole
ES =  initial_stuff_amounts(3);  %mole
P =  initial_stuff_amounts(4);   %mole

k1 = 2;   %1/(mole*second); rate of E and S binding 
k1r = 1;  %1/second; rate of ES disassociating
k2 = 1.5; %1/second; rate of product formation

dS_dt  = k1r*ES - k1*S*E;
dE_dt  = k1r*ES + k2*ES - k1*S*E;
dES_dt = k1*S*E - k1r*ES - k2*ES;
dP_dt  = k2*ES;

rates = [dS_dt; dE_dt; dES_dt; dP_dt];
end

%{
%description: explicit-eurler (failed attempt); the code is working but the plots may violate
some physics (e.g. negative amount of reagents) 
%parameters for rate:  
k1 = 2;   %1/(mole*second); rate of E and S binding 
k1r = 1;  %1/second; rate of ES disassociating
k2 = 1.5; %1/second; rate of product formation

%rate rules: 
syms S_curr;
syms E_curr;
syms ES_curr;
syms P_curr;
dS_dt  = k1r*ES_curr - k1*S_curr*E_curr;
dE_dt  = k1r*ES_curr + k2*ES_curr - k1*S_curr*E_curr;
dES_dt = k1*S_curr*E_curr - k1r*ES_curr - k2*ES_curr;
dP_dt  = k2*ES_curr;
       
%initial species amount: 
S_i =  8;   %mole
E_i =  4;   %mole
ES_i =  0;  %mole
P_i =  0;   %mole

h = 0.1; %deltal t
n = 100; %number of time points
t_intvl = linspace(0,10,101);

%% time evolution of the stuff amounts
S_t = zeros(1,n);
E_t = zeros(1,n);
ES_t = zeros(1,n);
P_t = zeros(1,n);                                     
S_t(1) = S_i;
E_t(1) = E_i;
ES_t(1) = ES_i;
P_t(1) = P_i;
for i =1:n
    S_t(:, i+1) = S_t(:,i) + h*subs(dS_dt,{S_curr, E_curr, ES_curr}, {S_t(:,i),E_t(:,i),ES_t(:,i)});
    E_t(:, i+1) = E_t(:,i) + h*subs(dE_dt,{S_curr, E_curr,ES_curr}, {S_t(:,i),E_t(:,i),ES_t(:,i)});
    ES_t(:, i+1) = ES_t(:,i) + h*subs(dES_dt,{S_curr, E_curr,ES_curr}, {S_t(:,i),E_t(:,i),ES_t(:,i)});
    P_t(:, i+1) = P_t(:,i) + h*subs(dP_dt,{ES_curr}, {ES_t(:,i)});
end
%}