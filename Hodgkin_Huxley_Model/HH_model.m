
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

%% Functions


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

