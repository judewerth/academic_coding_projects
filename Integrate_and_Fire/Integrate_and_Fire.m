%% Question 1
% Create Integrate and Fire neuron function
% Input: Time span vector, intial value (V0), Input Current
% Output: Voltage Vector
% Paarameters: Vrest, Vthresh, Vreset, Vpeak, Rm, Cm, Tref

%% b 
% Inputs
tspan = [-10, 100];
V0 = -60;


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