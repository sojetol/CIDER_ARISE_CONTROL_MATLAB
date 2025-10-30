% Playground


% t=0:0.001:1;
% B = 10;
% A = B* sin(2*pi*30*t);
% tau = 0.05;
% H = exp(-t);
% H1 = exp(-t/tau);
% y1 = 1 - H1;
% 
% y = A.*y1;
% 
% figure;
% plot(t,A)
% hold on
% plot(t,y)
% 



t = 0:0.001:1;
x = 10*sin(2*pi*60*t);          % input
tau = 0.002;
sys = tf(1, [tau 1]);           % H(s) = 1/(tau*s + 1)
y = lsim(sys, x, t);            % true LPF output

plot(t, x); hold on; plot(t, y); grid on;
legend('input','LPF output');
