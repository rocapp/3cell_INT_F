% Robert Capps, Dynamical Neuroscience assignment 1

clear all;
clf ;

% Define Parameters
I_Stim = 0.4;
R_m = 40; % Membrane resistance [MOhm]


t_f = [1.2 3.2 2.2]; % Time felt [s]
t_s = [3 2 1]; % Time sent[s]
alpha_1 =(t_f(1) - t_s(3)); % Time delay for cell 1
alpha_2 =(t_f(2) - t_s(1)); % Time delay for cell 2
alpha_3 =(t_f(3) - t_s(2)); % Time delay for cell 2
beta = 0.2;

% initial synaptic conductance
g_s = 0.017;


Erev=-15; % Inhibition
%Erev=15; % Excitation

dt = 1.; % Timestep [ms]


% Initialize Vectors
V_plot_vect=[];
S1_plot=[];
V_plot_vect_2=[];
S2_plot=[];
V_plot_vect_3=[];
S3_plot=[];





% Set initial Voltages
V_vect = -1; % First element of V. i.e. V at t=0
V_vect_2 = -2;
V_vect_3 = -3;

% Initial Conductance
S1=0;
S2=0;
S3=0;


t_end = 2500; % Total time [ms]

abs_ref = 10; % Absolute refractory period [ms]

ref1 = 0; % refractory period counter
ref2 = 0;
ref3 = 0;


V_th = 10; % Threshold voltage (spike) [mV]


    for t = 1:t_end %loop through values of t in steps of dt ms

    if ~ref1
    V_vect = (V_vect)-(V_vect/(R_m))+(I_Stim)-g_s*(V_vect-Erev)*(S3+S2);
     else
    ref1 = ref1 - 1;
    V_vect = -5;
    end

    if (V_vect > V_th) %cell spiked
        V_vect = 50; %set vector that will be plotted to show a spike here 
        ref1 = abs_ref;

    end
    
    if ~ref2
        V_vect_2 = (V_vect_2)-(V_vect_2/(R_m))+ I_Stim -g_s*(V_vect_2-Erev)*(S1+S3);
    else
        ref2 = ref2 - 1;
        V_vect_2 = -5;
    end

    % if statement below says what to do if voltage crosses threshold 
    if (V_vect_2 > V_th) %cell spiked
        V_vect_2 = 50; %set vector that will be plotted to show a spike here 
        ref2 = abs_ref;
    end

   

    
     if ~ref3
        V_vect_3 = (V_vect_3)-(V_vect_3/(R_m))+ I_Stim -g_s*(V_vect_3-Erev)*(S2+S1);
    else
        ref3 = ref3 - 1;
        V_vect_3 = -5;
    end

    %if statement below says what to do if voltage crosses threshold 
    if (V_vect_3 > V_th) %cell spiked
        V_vect_3 = 50; %set vector that will be plotted to show a spike here 
        ref3 = abs_ref;
    end

% loops through all time values, calculates synaptic conductance
for t= 1:t_end
  S1=S1+alpha_1*(1-S1)/(1+exp(-10*(V_vect-V_th)))-beta*S1;
  S2=S2+alpha_2*(1-S2)/(1+exp(-10*(V_vect_2-V_th)))-beta*S2;
  S3=S3+alpha_3*(1-S3)/(1+exp(-10*(V_vect_3-V_th)))-beta*S3;
end

% Initializes vectors
S3_plot = [S3_plot S3];
S2_plot = [S2_plot S2];
S1_plot=[S1_plot S1];

V_plot_vect = [V_plot_vect V_vect]; % With no spike, plot actual voltage V
V_plot_vect_2 = [V_plot_vect_2 V_vect_2];
V_plot_vect_3 = [V_plot_vect_3 V_vect_3];


end

% Plot Voltage

figure(1)
subplot(2,1,1); plot(V_plot_vect, 'Color', 'blue')
hold on;
title('Voltage vs. time')
subplot(2,1,1); plot(V_plot_vect_2, 'Color','green')
hold on;
subplot(2,1,1); plot(V_plot_vect_3, 'Color','r')
hold on;


xlabel('Time in ms');
ylabel('Voltage in mV');

whitebg([1 1 1])

% Plot Synaptic Current
subplot(2,1,2); plot(S1_plot, 'Color', 'blue')
hold on;
title('Synaptic Current vs. time')

subplot(2,1,2); plot(S2_plot, 'Color', 'green')
hold on;
subplot(2,1,2); plot(S3_plot, 'Color', 'red')
hold on;
whitebg([1 1 1])

% Find Phase lag

% Find the times when Cell 1 spikes
[pks,locs] = findpeaks(V_plot_vect,'MINPEAKHEIGHT',19);
SpikeTimes = locs


% Find the times when Cell 2 spikes
[pks_2,locs_2] = findpeaks(V_plot_vect_2,'MINPEAKHEIGHT',19);
SpikeTimes_2 = locs_2

% Find the times when Cell 3 spikes
[pks_3,locs_3] = findpeaks(V_plot_vect_3,'MINPEAKHEIGHT',19);
SpikeTimes_3 = locs_3

if length(SpikeTimes) > length(SpikeTimes_2)
    lastSpike = length(SpikeTimes_2);
else
    lastSpike = length(SpikeTimes);
end

if length(SpikeTimes) > length(SpikeTimes_3)
    lastSpike = length(SpikeTimes_3)
else
    lastSpike = length(SpikeTimes);
end

Period = [];
tau12 = [];
Tau12 = [];
tau13 = [];
Tau13 = [];

i = 1;

% Find period, tau
for i=1:lastSpike

if i < lastSpike
    i = i+1;
end

Period1 = (SpikeTimes(i) - SpikeTimes(i-1));
Period = [Period Period1];

tau12 = (SpikeTimes_2(i) - SpikeTimes(i-1));
tau13 = (SpikeTimes_3(i) - SpikeTimes(i-1));

Tau12 = [Tau12 tau12];
Tau13 = [Tau13 tau13];

end

% Find the transpose of the matrices
tTau12 = transpose(Tau12);
tTau13 = transpose(Tau13);
tPeriod = transpose(Period);

% Find Phase12 and Phase 13
Phase_12 = (tTau12/(tPeriod));

Phase_13 = (tTau13/(tPeriod));

% Modulo 1
Phi_13 = mod(Phase_13,1);
Phi_12 = mod(Phase_12,1);
    
% Length of the Phase matrices    

Spikes12 = length(Phase_12);
Spikes13 = length(Phase_13);

% Plot phases, Phi

figure (2)
subplot(3,1,1); plot(Phi_12, 'Color','blue')
xlabel('n (# of Spikes)');
ylabel('Phi_1_2');
axis([0 Spikes12 0 1]);
hold on;

subplot(3,1,2); plot(Phi_13, 'Color','red')
xlabel('n (# of Spikes)');
ylabel('Phi_1_3');
hold on;
axis([0 Spikes13 0 1]);
subplot(3,1,3); plot(Phi_12, Phi_13, 'Color','green')
xlabel('Phi_1_2');
ylabel('Phi_1_3');
hold on;
axis([0 1 0 1]);



