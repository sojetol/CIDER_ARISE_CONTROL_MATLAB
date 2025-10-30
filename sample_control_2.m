beta = 0.447;
tau = 1.946;
G = tf([beta],[tau 1])

[num,den]=pade(1,4);
s = tf('s');
G_delayed_pade = tf([beta],[tau 1])*tf(num,den)
G_delayed = tf([beta],[tau 1])*exp(-1*s)
G_delayed_less = tf([beta],[tau 1])*exp(-.25*s)

% s = tf("s");
% beta_d = 0.732
% tau_d = 4.063
% G_d = beta_d/(1+sqrt(tau_d*s))

figure
hold on
step(G)
step(G_delayed)
step(G_delayed_pade)
%%

opts = bodeoptions;
opts.PhaseWrapping = 'on';
figure
hold on
bodeplot(G,opts)
bodeplot(G_delayed_less,opts)
bodeplot(G_delayed,opts)
legend("No time delay","Time delay .25","Time delay 1")