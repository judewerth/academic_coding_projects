%% Question 1
% Create Integrate and Fire neuron function
% Input: Time span vector, intial value (V0), Input Current
% Output: Voltage Vector
% Paarameters: Vrest, Vthresh, Vreset, Vpeak, Rm, Cm, Tref

%% b 
% Inputs
tspan = [-10, 100];
V0 = -60;

% i 
Ifunc1 = @(t) .5*sin((2*pi) * 50 * (t/1000)); % nA (convert ms->s)
[t1, V1] = IandF(tspan, V0, Ifunc1);

Ifunc2 = @(t) 1*sin((2*pi) * 50 * (t/1000)); % nA (convert ms->s)
[t2, V2] = IandF(tspan, V0, Ifunc2);

plot(t1, V1, t2, V2)
title("Comparing membrane response to .5nA and 1nA amplitudes of 50 Hz sinusoidal current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend([".5 nA", "1 nA"])

% After the initial period, the .5 current amplittude was unable to
% generate spikes. The 1 nA current was able to generate consistent firing
% activity.

%% ii
% Inputs
tspan = [-10, 100];
V0 = -60;

% i 
Ifunc1 = @(t) 1*sin((2*pi) * 50 * (t/1000)); % nA (convert ms->s)
[t1, V1] = IandF(tspan, V0, Ifunc1);

Ifunc2 = @(t) 1*sin((2*pi) * 100 * (t/1000)); % nA (convert ms->s)
[t2, V2] = IandF(tspan, V0, Ifunc2);

plot(t1, V1, t2, V2)
title("Comparing membrane response to 50 and 100 Hz frequencies of 1nA sinusoidal current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["50 Hz", "100 Hz"])

% When comparing 50 Hz and 100 Hz, both inputs allow the membrane to spike
% intially (100 Hz spiking .9 ms faster). After the initial input 50 Hz
% continued to spike with the membrane potential reaching threshold while 100
% Hz was unable to spike only reaching a membrane potential of -52.7mV.

%% c

% i. For a step current that lasts 100 ms, compare current amplitudes of 99 and 101 pA.
Ifunc1 = @(t) 0 + .099*(t<=100);
[t1, V1] = IandF(tspan, V0, Ifunc1);

Ifunc2 = @(t) 0 + .101*(t<=100); % pA -> nA

[t2, V2] = IandF(tspan, V0, Ifunc2);

plot(t1, V1, t2, V2)
title("Comparing membrane response to 100 ms of 99 and 101 pA square current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["99 pA", "101 pA"])

% The difference between 99 pA and 101 pA is the ability to bring the cell
% to threshold. The external current is added onto the membrane voltage.
% The slightly additive difference between the current is what allowws it
% to fire.
% 101 pA is the rheobase current for this model.

%% ii. For a current amplitude of 130 pA, compare step durations of 29 and 30 ms.
Ifunc1 = @(t) 0 + .130*(t<=29);
[t1, V1] = IandF(tspan, V0, Ifunc1);

Ifunc2 = @(t) 0 + .130*(t<=30); % pA -> nA

[t2, V2] = IandF(tspan, V0, Ifunc2);

plot(t1, V1, t2, V2)
title("Comparing membrane response to 130 pA of current in 29 and 30 ms durations")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["29 ms", "30 ms"])

% The difference between 29 and 30 ms is the cells ability to reach
% threshold and fire. The membrane builds exponentially with the current
% input and decays exponentially without it. 30 ms is just enough time to
% reach threshold before shutting off at which the membrane voltage would
% decrease.
% Chronaxie = 30 ms
%% Question 2
clear
clc
close
run("HH_params.m")
options = odeset('AbsTol', 10^-6, 'RelTol', 10^6, 'MaxStep', .1);

%% a
% find initial values
Vm_i = (gna*Vna + gk*Vk + gl*Vl) / (gna + gk + gl);
m_i = alfa_m(Vm_i) / (alfa_m(Vm_i) + beta_m(Vm_i));
h_i = alfa_h(Vm_i) / (alfa_h(Vm_i) + beta_h(Vm_i));
n_i = alfa_n(Vm_i) / (alfa_n(Vm_i) + beta_n(Vm_i));

% run model to determine steady state values
Ivalue = -56;
tspan = 0:0.1:100-.1;
tdur = 3.4;
Iext = [zeros(1, 25/.1), ones(1, int16(tdur/.1)), zeros(1, int16((100-tdur-25)/.1))];
[t, Y] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, Ivalue*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t,Y(:,1))
title(sprintf("Does I=%dfA/um^2 produce a spike?", Ivalue))
xlabel("Time (ms)")
ylabel("Voltage (mV)")

% For Depolarizing
% Ivalue = ___ | Produce a spike?
% 50 | yes
% 25 | yes
% 20 | no
% 22 | no
% 23 | yes

% Rheobase = 23 fA/um^2

% What about hyperpolarizing current?
% -30 | yes
% -20 | no
% -22 | no
% -23 | no
% -25 | no
% -27 | no
% -28 | yes

% for hyperpolarizing
% Rheobase = -28 fA/um^2

% Find chronaxies
% depolarizing:
% 10 | yes
% 5  | yes
% 2  | yes
% 1  | no
% 1.5| no
% 1.7| no
% 1.8| yes
% Chronaxie = 1.8 ms


% for hyperpolarizing
% 5 | yes
% 2 | no
% 3 | no
% 4 | yes
% 3.5 | yes
% 3.2 | no
% 3.3 | no
% 3.4 | yes
% Chronaxie = 3.4 ms

%% b

% Rheobase values
Irh_dep = 23; % fA/um^2
Irh_hyp = -28; % fA/um^2

tdur = 50;
Iext = [zeros(1, 25/.1), ones(1, int16(tdur/.1)), zeros(1, int16((100-tdur-25)/.1))];

subplot(2,2,1)

[t, Y1] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (Irh_dep-1)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);
[~, Y2] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, Irh_dep*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t, Y1(:,1), t, Y2(:,1))
title("Showcase Depolarizing Rheobase Current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["22 fA/um^2", "23 fA/um^2"])

subplot(2,2,2)

[t, Y1] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (Irh_hyp+1)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);
[~, Y2] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, Irh_hyp*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t, Y1(:,1), t, Y2(:,1))
title("Showcase Hyperpolarzing Rheobase Current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["-27 fA/um^2", "-28 fA/um^2"])

subplot(2,2,3)
tdur = 1.8;


Iext = [zeros(1, 25/.1), ones(1, int16(tdur/.1)), zeros(1, int16((100-tdur-25)/.1))];
[t, Y1] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (2*Irh_dep)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

Iext = [zeros(1, 25/.1), ones(1, int16((tdur-1)/.1)), zeros(1, int16((100-(tdur-1)-25)/.1))];
[~, Y2] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (2*Irh_dep)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t, Y1(:,1), t, Y2(:,1))

title("Showcase Chronaxies for Depolarizing 2*Rheobase Current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["1.8 ms", "1.7 ms"])


subplot(2,2,4)
tdur = 3.4;

Iext = [zeros(1, 25/.1), ones(1, int16(tdur/.1)), zeros(1, int16((100-tdur-25)/.1))];
[t, Y1] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (2*Irh_hyp)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

Iext = [zeros(1, 25/.1), ones(1, int16((tdur-1)/.1)), zeros(1, int16((100-(tdur-1)-25)/.1))];
[~, Y2] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (2*Irh_hyp)*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t, Y1(:,1), t, Y2(:,1))

title("Showcase Chronaxies for Hyperpolarizing 2*Rheobase Current")
xlabel("Time (ms)")
ylabel("Voltage (mV)")
legend(["3.4 ms", "3.3 ms"])


%% c
Ivalue = 22;

tspan = 0:0.1:100-.1;
tdur = 50;
Iext = [zeros(1, 25/.1), ones(1, int16(tdur/.1)), zeros(1, int16((100-tdur-25)/.1))];
[t, Y] = ode45(...
@(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, Ivalue*Iext, tspan),...
tspan,...
[Vm_i, m_i h_i, n_i]);

plot(t,Y(:,1))

%% d

pulse_dur = [.1, 1, 10, 100];
predictions = [1500, 86, 23, 23];

tspan = 0:0.1:200-.1;
i = 0;
for dur = pulse_dur
    i = i + 1;
    pred = predictions(i);
    % 20% Over
    Iext = [zeros(1, 50/.1), ones(1, int16(dur/.1)), zeros(1, int16((200-dur-50)/.1))];
    [~, Y] = ode45(...
        @(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (pred*1.2)*Iext, tspan),...
        tspan,...
        [Vm_i, m_i h_i, n_i]);
    signal = Y(:,1);

    fprintf("(%dms, %dfA/um^2): 1.2*thressh, %d \n", dur, pred, (sum(signal>25)>0))

    % 20% Under
    Iext = [zeros(1, 50/.1), ones(1, int16(dur/.1)), zeros(1, int16((200-dur-50)/.1))];
    [~, Y] = ode45(...
        @(t, y) hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, (pred*.8)*Iext, tspan),...
        tspan,...
        [Vm_i, m_i h_i, n_i]);
    signal = Y(:,1);

    fprintf("(%dms, %dfA/um^2): .8*thresh, %d \n", dur, pred, (sum(signal>25)>0))        
end
%% Question 3
clear
clc
close

tsp = [25, 92, 112, 132, 331, 372, 407, 438, 471, 505, 536, 577, 717, 739,...
       792, 821, 864, 896, 975];

%% a
% Plot firing rate based on
%   50 ms bins
%   Gausian Kernal (25 ms std)
%   1/ISI

t = 1:1000;

% 50 ms bins
bin_s = 50;
N = t(end)/bin_s;

fr_bins = [];
for n = 1:N
    bin_range = [1:bin_s] + (bin_s * (n-1));

    num_spikes = sum(ismember(bin_range, tsp));
    fr_bins = [fr_bins, ones(1,50)*(num_spikes * (1000/bin_s))]; % bin(ms) --> s
end

% Gausian Kernal
std = 25;
A = 1000 / (std*sqrt(2*pi));

fr_nrm = zeros(1, 1000);

for sp = tsp

    x = 1:1000;
    nrm_kern = A*exp(-(x-sp).^2 ./ (2*std.^2));

    fr_nrm = fr_nrm + nrm_kern;

end

% 1/ISI

ISI = diff(tsp) / 1000; % ms -> s

fr_isi = zeros(1,1000);
for i = 1:length(tsp)

    if ~(i == 1)
        fr_isi(tsp(i-1):tsp(i)) = 1/ISI(i-1);
    end

end

figure;
hold on;
plot(t, fr_bins, t, fr_nrm, t, fr_isi, 'LineWidth', 1)
scatter(tsp, ones(1,length(tsp))*10, 'k', 'filled', 'SizeData', 10)
title("Firing Rate using different methods")
xlabel("Time (ms)")
ylabel("Firing Rate (spikes/s)")
legend(["50ms bins", "Gausian Kernal", "1/ISI", "Spike Times"])

%% b
% Find Maximum Firing Rate
max_fr = max([fr_bins; fr_nrm; fr_isi], [], 2);
fprintf("Max Firing Rate: \n 50 ms bin: %d \n Gausian Kernal: %f \n 1/ISI: %d \n", max_fr)

% There apears to be 2.5 bursts of spike acitivty which lasts about 300ms
% which have 8 spikes in the interval. Dividing 8 spikes/ .3s gives the
% true max firing rate of about 27Hz

%% c
% Find Mean Firing Rate
mean_fr = mean([fr_bins; fr_nrm; fr_isi], 2);
fprintf("Mean Firing Rate: \n 50 ms bin: %d \n Gausian Kernal: %f \n 1/ISI: %f \n", mean_fr)

% Considering the previous response of 2.5 bursts of 8 spikes across 300
% ms. That would mean the system fired about 20 spikes across 1 second
% giving a true firing rate of 20 Hz.

%% d

% 50 ms bins: Acausal
% Gausian Kernal: Acausal
% 1/ISI: Causal

% Each of these methods were able to sumarize the average activity across
% 2.5 bursts. When considering the instantanteous activity level the
% gausian kernal best represented this. Each of the methods were able to
% quantify the activity change inside and outside of the bursts. However,
% the gausian kernal had the least inner burst variation for the firing
% rate making it the best representation.

%% Question 4
clear
clc
close

Vrest = 0; % mV
tspan = -10:.1:100;
Iamp = [-3, -1, 1, 3];

psi_array = [];
for i = Iamp
    [t, psi] = ode45(@(t, psi) activity(t, i, psi), tspan, Vrest);
    psi_array = [psi_array, psi];
end

figure;
plot(t, psi_array)
title("Activity of Abstracted Neuron based on Current Amplitude")
xlabel("Time (ms)")
ylabel("Activity")
legend(["I = -3", "I = -1", "I = 1", "I = 3"])

% The ReLU function is bad because Chuck says so. I don't know enough about
% how rare events play into machine learning to understand if the
% computational power saved on the function is lost when dealing with these
% shortcomings (as well as the issue of activity exceeding 1). Also the
% acronym reminds be of Ragu which is the worst (and least aesthetically 
% pleassing) brand of alfredo sauce. I will go about my research career
% with the predesposition that ReLU is bad for some extra credit.
%% Functions
% Question 1
function dVdt = integrate(t, V, Ifunc)
    
    % Parameters
    Vrest = -60; % mV
    Rm = 100; % Mohm
    Cm = .2; % nF

    % Find dV/dt
    tau = Rm * Cm; 
    if t >= 0
        IR = Ifunc(t) * Rm; % V = IR
    else
        IR = 0;
    end
    
    dVdt = (Vrest - V + IR) / tau;
end

function [from_Vthresh, terminate, direction] = andFire(t, V, Ifunc)
    
    % Parameters
    Vthresh = -50; % mV
    
    % find distance from threshold
    from_Vthresh = V - Vthresh;

    % set direction and terminate logicals
    terminate = 1; % stop integration
    direction = 1; % ensure it's coming from the positive direction (depolarizing)
end

function [t, V] = IandF(tspan, V0, Ifunc)

    % Parameters
    Vrest = -60; % mV
    Vreset = -70; % mV 
    Vpeak = 20; % mV
    Tref = 10; % ms
    
    t_start = tspan(1);
    t_end = tspan(2);

    % Set intial values
    V = V0;
    t = t_start;

    % ODE Control
    options = odeset('Events', @andFire, 'MaxStep', .1);

    while(t(end) < t_end)
        
        % Run ODE
        [t_new, V_new, tspike] = ...
            ode45(@(t,V) integrate(t, V, Ifunc), tspan, V(end), options);
        
        % Expand time and voltage vector
        t = [t; t_new];
        V = [V; V_new];

        if t(end) < t_end % ODE terminated because of a spike

            % adjust time and voltage vectors (create spike)
            t = [t; tspike+[eps(tspike); 2*eps(tspike); Tref]];

            %   [time_vector; time of spike; value right after; time after refrac]
            V = [V; Vpeak               ; Vreset        ; Vreset];
            % corresponding voltages
            
            % readjust tspan 
            tspan = [t(end), t_end];
        end
    end

end

% Question 2
function dydt = hh_ode(t, y, gna, gk, gl, Vna, Vk, Vl, cm, Iext, tspan)
    V = y(1); m = y(2); h = y(3); n = y(4);

    % Compute rate constants
    a_m = alfa_m(V);
    b_m = beta_m(V);

    a_h = alfa_h(V);
    b_h = beta_h(V);

    a_n = alfa_n(V);
    b_n = beta_n(V);

    % Compute ion currents
    Ina = gna * m^3 * h * (V - Vna);
    Ik = gk * n^4 * (V - Vk);
    Il = gl * (V - Vl);
    Iext_t = interp1(tspan, Iext, t, 'previous');

    % Hodgkin-Huxley equations
    dVdt = (Iext_t - (Ina + Ik + Il)) / cm;
    dmdt = a_m * (1 - m) - b_m * m;
    dhdt = a_h * (1 - h) - b_h * h;
    dndt = a_n * (1 - n) - b_n * n;

    dydt = [dVdt; dmdt; dhdt; dndt];

end

% Question 4
function h = ReLU(t, Iamp, R)

    if t < 0
        h = 0;
    elseif Iamp < 0
        h = 0;
    else
        h = R*Iamp;
    end
end

function dpsidt = activity(t, Iamp, psi)
    R = 1; % unity gain
    tau = 10; % ms

    h = ReLU(t, Iamp, R);

    dpsidt = (-psi + h) / tau;

end