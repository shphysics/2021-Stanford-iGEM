clear all; close all; clc; 
%{
curr version: v1 
1. standard yeast growth curve 
2. no plasmid-burden on yeast colony
3. assume existing CYP will not affect future CYP prudction rate
4. assume CYP are not degrading 
%}

%% basic constants
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


%% plotting
figure(1);
hold on; 
plot(t_intvl, subs(CYP,t_intvl));
plot(t_intvl, subs(CYP1,t_intvl));
plot(t_intvl, subs(CYP2,t_intvl));
xlabel('time (hours)');
ylabel('amounts (moles)');
legend('CYP amounts B2=0.3,B3=2', 'CYP amounts B2=0.9,B3=6', 'CYP amounts B2=0.1,B3=12');
title('CYP production curve');
hold off; 