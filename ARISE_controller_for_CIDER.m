% Explicit feedback for climate modeling
% Copyright (C) 2020  Ben Kravitz
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.

% Sample control parameters file
%
% Written by Ben Kravitz (bkravitz@iu.edu or ben.kravitz.work@gmail.com)
% Last updated 11 July 2019
%
% This script provides information about the feedback algorithm.  All of this
% is user-customizable.  The other parts of the script will give you outvals,
% lats, lons, and times.  The output of this script should be a list called
% nlvals, which consists of pairs.  The first item in the pair is the name
% of the namelist value.  The second item is the value itself as a string.

% Translated to MATLAB by Sandia AI and some fiddling around with it by
% Jared Farley to make it so the parameters are an input into the function

function [q,log_array] = ARISE_controller_for_CIDER(outvals,log_array,controller_params,preset_noise)
% USER-SPECIFIED CONTROL PARAMETERS
% refvals = [288.64, 0.8767, -5.89]; % updated to be average over years 2010-2029


refvals = controller_params.refvals; % T0 T1 T2 Targets
kivals = controller_params.kivals; % Integral Gains
kpvals = controller_params.kpvals; % Proportional gains
summing_index = controller_params.months_to_average-1; % How long to average over to determine "current" values
sampling_period = controller_params.sampling_period; % how long since last time controller was called

% new_refvals = [288.14, 0.8497, -5.90]; % new target values: averages during the years 2000-2019
% new_refvals = [287.64, 0.7348, -5.97]; % new target values: averages during the years 1984-2003
% refvals = [288.21, 0.594, -6.006]; % new version of the model (GLENS values)
% kivals = [0.0183, 0.0753, 0.3120];
% kpvals = [0.0183, 0.0753, 0.3120];


firstyear = 2035;
baseyear = 2030;
% x_ramp = 5.0; % defines a range of years over which the feedback is ramped up
x_ramp = controller_params.x_ramp; % defines a range of years over which the feedback is ramped up
size_outvals = size(outvals);

lats = linspace(-90,90,size_outvals(2));
lons = linspace(0,360,size_outvals(1));

% USER SPECIFIED CALCULATIONS

% logheader = {'Timestamp', 'dT0', 'sum(dT0)', 'dT1', 'sum(dT1)', 'dT2', 'sum(dT2)', 'L0', 'L1N', 'L1S', 'L2', '30S(Tg)', '15S(Tg)', '15N(Tg)', '30N(Tg)'};

firsttime = 0;
if isempty(log_array)
    firsttime = 1;
else
    loglines = log_array;
end


% w = makeweights(lats', lons');
T0 = mean(globalMean(outvals(:,:,end-summing_index:end))+preset_noise(end-summing_index:end,1));
T1 = mean(calculateT1(outvals(:,:,end-summing_index:end))+preset_noise(end-summing_index:end,2));
T2 = mean(calculateT2(outvals(:,:,end-summing_index:end))+preset_noise(end-summing_index:end,3));

de = [T0 - refvals(1), T1 - refvals(2), T2 - refvals(3)]; % error terms

if firsttime == 1
    timestamp = firstyear;
    sumde = de/12*sampling_period;
    sumdt2 = de(3)/12*sampling_period;
else
    timestamp = loglines(end, 1) + sampling_period/12;
    sumdt0 = loglines(end, 3) + (T0 - refvals(1))/12*sampling_period;
    sumdt1 = loglines(end, 5) + (T1 - refvals(2))/12*sampling_period;
    sumdt2 = loglines(end, 7) + (T2 - refvals(3))/12*sampling_period;
    sumde = [sumdt0, sumdt1, sumdt2];
end

dt = timestamp - baseyear;
dt2 = timestamp - firstyear;

sens=4.1;
sens=controller_params.sens;
change=controller_params.ff_rate*dt+(288.5029-refvals(1)); 

l0hat=change/sens/1.40;
% l0hat = 0;
% Updated based on feedback simulation
% l0hat = 0.00347 * dt;
l1hat = 0.000 * dt;
l2hat = 0.00 * dt;

ramp_up = 1.0;
if (dt2 < x_ramp)
    ramp_up = dt2 / x_ramp;
end

% Feedback
l0kp1 = (kpvals(1) * de(1) + kivals(1) * sumde(1)) * ramp_up;
l1kp1 = (kpvals(2) * de(2) + kivals(2) * sumde(2) - 0.5 * l0kp1) * ramp_up;
l2kp1 = (kpvals(3) * de(3) + kivals(3) * sumde(3) - l0kp1) * ramp_up;

% All of the feeds
l0step4 = l0kp1 + l0hat;
% [l0kp1 l0hat]
l1step4 = l1kp1 + l1hat;
l2step4 = l2kp1 + l2hat;

l0 = max(l0step4, 0);
l1n = min(max(l1step4, 0), l0);
l1s = min(max(-l1step4, 0), l0);
l2 = min(max(l2step4, 0), l0 - l1s - l1n);
ell = [l0; l1n; l1s; l2];

% Preventing integrator wind-up
if (l2 == (l0 - l1s - l1n))
    sumdt2 = sumdt2 - (T2 - refvals(3));
    sumde(3) = sumdt2;
end



M = [0, 30, 30, 0; 
     0, 0, 45, 20; 
     20, 45, 0, 0; 
     40, 0, 0, 40];
M = controller_params.M;
F = [1, 1, 1, 1; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];

% M = [0 35.5 35.7 0;
%      0 10.8 45.5 18.8;
%      17.3 46.9 10.5 0;
%      39.1 6.6 4.8 39.7];

q = M' / F * ell;

for k = 1:length(q)
    q(k) = max(q(k), 0);
end

% newline = {num2str(timestamp), num2str(de(1)), num2str(sumde(1)), num2str(de(2)), num2str(sumde(2)), num2str(de(3)), num2str(sumde(3)), num2str(l0), num2str(l1n), num2str(l1s), num2str(l2), num2str(q(1)), num2str(q(2)), num2str(q(3)), num2str(q(4))};
newline = [timestamp, de(1), sumde(1), de(2), sumde(2), de(3), sumde(3), l0, l1n, l1s, l2, q(1), q(2), q(3), q(4)];

log_array = [log_array;newline];

% if firsttime == 1
%     linestowrite = [logheader; newline];
% else
%     linestowrite = loglines;
%     linestowrite(end + 1, :) = newline;
% end
% 
% writelog(fullfile(maindir, logfilename), linestowrite);

% USER SPECIFIED OUTPUT
end