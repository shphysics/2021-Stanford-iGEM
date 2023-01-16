%Overview of growth_curve 
%{
the growth model should help to construct growth curve for both the un-engineered type and engineered type yeasts 
this allows us to see if incorporating the two plasmids into the yeast will
affect its growth ability

based on paper: http://www.math.chalmers.se/Stat/Research/Preprints/Doctoral/2005/3.pdf
%}
clear all; close all; clc; 

%% setting up time interval (in hours)
t_intvl = linspace(0,60,61);

%% setting up parameters (parameters values are based on paper Pylvanainen (2005))
syms t; % place holder for time 

%parameters with physical meanings (expressed as tunable lower-level parameters)

%{
B0 = Az; %the asymptote, the maximum value of the growth reached (on the log
scale)

B0*B2*((B3)^(B3/(1-B3)))= growth_rate; %max rel. pop. growth rate
 
(B0(1-B1)^(1/(1-B3)) - B0*B3^(1/(1-B3)) + growth_rate*(log(B1/(1-B3))/B2))/growth_rate= lag; %time for the yeast to adapt to env. before exp. growth

(log(B1/(1-B3)))/B2= time_I; %inflection time point 
%}

% lower-level parameters (by obtaining the physical parameters above from
% experiment, we can solve for the lower-level parameters 
B0 = 4.5; 
B1 = -50;
B2 = 0.3;
B3 = 3; % 2/3: von Bertalanffy function; 2: logistic model 
D = -3; 

B01 = 6;
B02 = 2;

B11 = -30;
B12 = -70;

B21 = 0.9;
B22 = 0.1;

B31 = 6;
B32 = 12;
%% champman-Richards model (B0,B2 > 0, B3 > 1, and B1 < 1 - B3 or B0,B2 > 0, 0 < B3 < 1, and 1 - B3 < B1 < 1)
gt = B0*(1 - (B1*exp((-B2)*t)))^(1/(1 - B3)) + D; %log(Nt);
gt1 = B0*(1 - (B1*exp((-B21)*t)))^(1/(1 - B31)) + D; %log(Nt);
gt2 = B0*(1 - (B1*exp((-B22)*t)))^(1/(1 - B32)) + D; %log(Nt);
%N_t = 10^(gt); %the Nt; 

val_gt = subs(gt, t_intvl);

%% plotting C-R model with 
figure(1); 
hold on; 
plot(t_intvl, subs(gt, t_intvl));
plot(t_intvl, subs(gt1, t_intvl));
plot(t_intvl, subs(gt2, t_intvl));
xlabel('time (hours)');
ylabel('log of population size (log(N_t))');
legend('B3=3', 'B3=6', 'B3=12');
title('Saccharomyces cerevisiae predicted growth curve');
hold off; 
