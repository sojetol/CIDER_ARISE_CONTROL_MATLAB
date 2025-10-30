%% Variables to change
% System is a block with constant and oscilatory forces acting on it
% We want to make it stay in a spot
my_title = "PD Control with different sampling rates, k_p = 0.25, w_{osc} = 5";
% my_title = "PD Control with different sampling rates, k_p = 0.25, no oscillatory force";

f = 1; % Constant force on the block
f0 = 0; % Constant force on the block
ang_velo = 5; % Angular velo of the oscilatory force on the block
freq_noise = "11*pi";
amp = 50; % amplitude of oscilatory force on block (ish)
%amp = 0; % amplitude of oscilatory force on block (ish)
x_target = -10; % Targeted position
kp = .25; % Proportional gain
kd = .5; % Derivative gain 
holder = 0; % Placeholder for annual or seasonal position information
dholder = 0;  % Placeholder for annual or seasonal velocity information

vars = [f;ang_velo;kp;holder;kd;dholder;amp;x_target];
vars0 = [f0;ang_velo;kp;holder;kd;dholder;amp;x_target];


end_time = 50; % Number of 'years'
points_in_year = 120; % Data points in year (ensure div by 4)
%% Code
tseries = linspace(0,end_time,end_time*points_in_year+1);
% plant = sin(pi*t)+f*t.^2/10;%+ sin(2*pi*t+pi/4);
[t_plant,plant] = ode45(@(t,y) dydt_uncontrolled(t,y,vars),tseries,[0;amp*ang_velo]);
[t_plant,plant0] = ode45(@(t,y) dydt_uncontrolled(t,y,vars0),tseries,[0;amp*ang_velo]);
[t_plant,controlled_plant] = ode45(@(t,y) dydt(t,y,vars),tseries,[0;amp*ang_velo]);
t_plant2 = zeros(points_in_year*end_time,1);
t_plant3 = zeros(points_in_year*end_time,1);
what_controller_sees = [];
what_controller_sees_2 = [];
what_controller_sees_3 = [];
what_controller_sees_4 = [];
controlled_plant_sparse = zeros(end_time*points_in_year,2);
controlled_plant_sparse_zohseas = zeros(end_time*points_in_year,2);
controlled_plant_sparse_seas = zeros(end_time*points_in_year,2);
for i = 0:(end_time-1)
    if i == 0
        [t_plant2(points_in_year*i+1:points_in_year*i+points_in_year),controlled_plant_sparse(points_in_year*i+1:points_in_year*i+points_in_year,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year*i+1:points_in_year*i+points_in_year),[0;amp*ang_velo]);
        
        vars(4) = mean(controlled_plant_sparse(points_in_year*i+1:points_in_year*i+points_in_year,1));
        what_controller_sees = [what_controller_sees;vars(4)];
    else
        [t_plant2(points_in_year*i:points_in_year*i+points_in_year),controlled_plant_sparse(points_in_year*i:points_in_year*i+points_in_year,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year*i:points_in_year*i+points_in_year),controlled_plant_sparse(points_in_year*i,:)');
        vars(4) = mean(controlled_plant_sparse(points_in_year*i+1:points_in_year*i+points_in_year,1));
        vars(6) = (vars(4)-mean(controlled_plant_sparse(points_in_year*(i-1)+1:points_in_year*(i-1)+points_in_year,1)))/1;
        what_controller_sees = [what_controller_sees;vars(4)];
    end
end
vars(4) = 0;
vars(6) = 0;
for i = 0:(end_time-1)
    if i == 0
        [t_plant2(points_in_year*i+1:points_in_year*i+points_in_year),controlled_plant_sparse_zohseas(points_in_year*i+1:points_in_year*i+points_in_year,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year*i+1:points_in_year*i+points_in_year),[0;amp*ang_velo]);

        vars(4) = mean(controlled_plant_sparse_zohseas(points_in_year*i+1:points_in_year*i+points_in_year,1));
        what_controller_sees_4 = [what_controller_sees_4;vars(4)];
    else
        [t_plant2(points_in_year*i:points_in_year*i+points_in_year),controlled_plant_sparse_zohseas(points_in_year*i:points_in_year*i+points_in_year,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year*i:points_in_year*i+points_in_year),controlled_plant_sparse_zohseas(points_in_year*i,:)');
        vars(4) = mean(controlled_plant_sparse_zohseas(points_in_year*i+1+3/4*points_in_year:points_in_year*i+points_in_year,1));
        two_seasons_ago = mean(controlled_plant_sparse_zohseas(points_in_year*i+1+2/4*points_in_year:points_in_year*i+3/4*points_in_year,1));
        vars(6) = (vars(4)-two_seasons_ago)/(1/4);
        what_controller_sees_4 = [what_controller_sees_4;vars(4)];
    end
end
vars(4) = 0;
vars(6) = 0;
for i = 0:(end_time*4-1)
    if i == 0
        [t_plant2(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4),[0;amp*ang_velo]);

        vars(4) = mean(controlled_plant_sparse_seas(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,1));
        what_controller_sees_2 = [what_controller_sees_2;vars(4)];
    else
        [t_plant2(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas(points_in_year/4*i:points_in_year/4*i+points_in_year/4,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas(points_in_year/4*i,:)');
        vars(4) = mean(controlled_plant_sparse_seas(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,1));
        vars(6) = (vars(4)-mean(controlled_plant_sparse_seas(points_in_year/4*(i-1)+1:points_in_year/4*(i-1)+points_in_year/4,1)))/(1/4);
        what_controller_sees_2 = [what_controller_sees_2;vars(4)];
    end
end
vars(4) = 0;
vars(6) = 0;
for i = 0:(end_time*4-1)
     if i == 0
        [t_plant3(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas_filt(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4),[0;amp*ang_velo]);

        vars(4) = mean(controlled_plant_sparse_seas_filt(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,1));
        what_controller_sees_3 = [what_controller_sees_3;vars(4)];
     elseif i<4
        [t_plant3(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas_filt(points_in_year/4*i:points_in_year/4*i+points_in_year/4,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas_filt(points_in_year/4*i,:)');

        vars(4) = mean(controlled_plant_sparse_seas_filt(points_in_year/4*i+1:points_in_year/4*i+points_in_year/4,1));
        what_controller_sees_3 = [what_controller_sees_3;vars(4)];
    else
        [t_plant3(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas_filt(points_in_year/4*i:points_in_year/4*i+points_in_year/4,:)] = ode45(@(t,y) dydt_sparse(t,y,vars),tseries(points_in_year/4*i:points_in_year/4*i+points_in_year/4),controlled_plant_sparse_seas_filt(points_in_year/4*i,:)');
        vars(4) = mean(controlled_plant_sparse_seas_filt(points_in_year/4*(i-3)+1:points_in_year/4*i+points_in_year/4,1));
        vars(6) = (vars(4)-mean(controlled_plant_sparse_seas_filt(points_in_year/4*(i-4)+1:points_in_year/4*(i-1)+points_in_year/4,1)))/(1/4);
        what_controller_sees_3 = [what_controller_sees_3;vars(4)];
    end
end
%%
figure
hold on
box on 
grid on
line([5 end_time],[x_target x_target],LineWidth=2,LineStyle="--",Color='g')
plot(t_plant,plant(:,1),LineWidth=2,Color='k')
plot(t_plant2,controlled_plant_sparse(:,1),LineWidth=2,Color='r')
plot(t_plant2,controlled_plant_sparse_seas(:,1),LineWidth=2,Color='b')
plot(t_plant2,controlled_plant_sparse_seas_filt(:,1),LineWidth=2,Color='c')
% plot(t_plant,controlled_plant(:,1),LineWidth=2)
% plot(t_plant,plant(:,1),LineWidth=2,Color='r')
ylim([-60 200])

sparse_x = reshape((what_controller_sees.*ones(end_time,points_in_year))',points_in_year*(end_time),[]);
sparse_x_seas = reshape((what_controller_sees_2.*ones(end_time*4,points_in_year/4))',points_in_year/4*(end_time*4),[]);
sparse_x_seas_filt = reshape((what_controller_sees_3.*ones(end_time*4,points_in_year/4))',points_in_year/4*(end_time*4),[]);
plot(t_plant2(2*points_in_year+1:end),sparse_x(points_in_year+1:end-points_in_year),LineWidth=2,LineStyle="--",Color='r')
plot(t_plant2(2*points_in_year/4+1:end),sparse_x_seas(points_in_year/4+1:end-points_in_year/4),LineWidth=2,LineStyle="--",Color='b')
plot(t_plant2(2*points_in_year/4+1:end),sparse_x_seas_filt(points_in_year/4+1:end-points_in_year/4),LineWidth=2,LineStyle="--",Color='c')
title(my_title)
xlabel("Time (year)")
ylabel("Position")
legend("Target","Uncontrolled","Last-Year-Mean Control","Seasonal Control","Seasonal Control w/ LYM Filter","What LYM Controller Sees","What Seasonal Controller Sees","What Seasonal Controller w/ LYM Filter Sees","Location","nw")

%%
figure
hold on
box on 
grid on
line([5 end_time],[x_target x_target],LineWidth=2,LineStyle="--",Color='g')
plot(t_plant,plant(:,1),LineWidth=2,Color='k')
plot(t_plant2,controlled_plant_sparse(:,1),LineWidth=2,Color='r')
plot(t_plant2,controlled_plant_sparse_zohseas(:,1),LineWidth=4,Color='m')
plot(t_plant2,controlled_plant_sparse_seas(:,1),LineWidth=2,Color='b')
plot(t_plant2,controlled_plant_sparse_seas_filt(:,1),LineWidth=2,Color='c')
plot(t_plant,controlled_plant(:,1),"Color",'#00cc88',LineWidth=2)
% plot(t_plant,plant(:,1),LineWidth=2,Color='r')

sparse_x = reshape((what_controller_sees.*ones(end_time,points_in_year))',points_in_year*(end_time),[]);
sparse_x_seas = reshape((what_controller_sees_2.*ones(end_time*4,points_in_year/4))',points_in_year/4*(end_time*4),[]);
sparse_x_seas_filt = reshape((what_controller_sees_3.*ones(end_time*4,points_in_year/4))',points_in_year/4*(end_time*4),[]);
% plot(t_plant2(2*points_in_year+1:end),sparse_x(1:end-points_in_year),LineWidth=2)
% plot(t_plant2(2*points_in_year/4+1:end),sparse_x_seas(1:end-points_in_year/4),LineWidth=2)
title(my_title)
xlabel("Time (year)")
ylabel("Position")
ylim([-60 200])

legend("Target","Uncontrolled","Last-Year-Mean Control","Last-Season-Held-For-a-Year Control","Seasonal Control","Seasonal Control w/ LYM Filter","Instantaneous Control","Location","nw")

%%
figure
hold on
box on 
grid on
line([5 end_time],[x_target x_target],LineWidth=2,LineStyle="--",Color='g')
plot(t_plant,plant(:,1)-plant0(:,1),LineWidth=2,Color='k')
plot(t_plant2,controlled_plant_sparse(:,1)-plant0(1:end-1,1),LineWidth=2,Color='r')
plot(t_plant2,controlled_plant_sparse_zohseas(:,1)-plant0(1:end-1,1),LineWidth=4,Color='m')
plot(t_plant2,controlled_plant_sparse_seas(:,1)-plant0(1:end-1,1),LineWidth=2,Color='b')
plot(t_plant2,controlled_plant_sparse_seas_filt(:,1)-plant0(1:end-1,1),LineWidth=2,Color='c')
plot(t_plant,controlled_plant(:,1)-plant0(:,1),"Color",'#00cc88',LineWidth=2)
% plot(t_plant,plant(:,1),LineWidth=2,Color='r')

sparse_x = reshape((what_controller_sees.*ones(end_time,points_in_year))',points_in_year*(end_time),[]);
sparse_x_seas = reshape((what_controller_sees_2.*ones(end_time*4,points_in_year/4))',points_in_year/4*(end_time*4),[]);
% plot(t_plant2(2*points_in_year+1:end),sparse_x(1:end-points_in_year),LineWidth=2)
% plot(t_plant2(2*points_in_year/4+1:end),sparse_x_seas(1:end-points_in_year/4),LineWidth=2)
title(my_title)
xlabel("Time (year)")
ylabel("Position")
ylim([-60 200])
title(my_title+", Noise-removed")
legend("Target","Uncontrolled","Last-Year-Mean Control","Last-Season-Held-For-a-Year Control","Seasonal Control","Seasonal Control w/ LYM Filter","Instantaneous Control","Location","nw")
%%
% plant_no_noise = plant;
% controlled_plant_sparse_no_noise = controlled_plant_sparse;
% controlled_plant_sparse_zohseas_no_noise = controlled_plant_sparse_zohseas;
% controlled_plant_sparse_seas_no_noise = controlled_plant_sparse_seas;
% controlled_plant_sparse_seas_filt_no_noise = controlled_plant_sparse_seas_filt;
% controlled_plant_no_noise = controlled_plant;
%%
figure
hold on
box on 
grid on
% line([5 end_time],[x_target x_target],LineWidth=2,LineStyle="--",Color='g')
plot(t_plant,plant(:,1)-plant0(:,1)-plant_no_noise(:,1),LineWidth=2,Color='k')
plot(t_plant,controlled_plant(:,1)-plant0(:,1)-controlled_plant_no_noise(:,1),"Color",'#00cc88',LineWidth=2)
plot(t_plant2,controlled_plant_sparse_seas(:,1)-plant0(1:end-1,1)-controlled_plant_sparse_seas_no_noise(:,1),LineWidth=2,Color='b')
plot(t_plant2,controlled_plant_sparse_seas_filt(:,1)-plant0(1:end-1,1)-controlled_plant_sparse_seas_filt_no_noise(:,1),LineWidth=2,Color='c')
plot(t_plant2,controlled_plant_sparse(:,1)-plant0(1:end-1,1)-controlled_plant_sparse_no_noise(:,1),LineWidth=2,Color='r')
plot(t_plant2,controlled_plant_sparse_zohseas(:,1)-plant0(1:end-1,1)-controlled_plant_sparse_zohseas_no_noise(:,1),LineWidth=2,Color='m')
title(my_title)
xlabel("Time (year)")
ylabel("Position")
title("Effect of Noise On Behavior, w_{noise} = "+freq_noise)
legend("0","Uncontrolled (numerical drift)","Last-Year-Mean Control","Seasonal Control","Instantaneous Control","Location","nw")
legend("Uncontrolled (numerical drift)","Instantaneous Control","Seasonal Control (8pi)","Seasonal Control w/ LYM Filter (8pi)","Last-Year-Mean Control (2pi)","Last-Season-Held-For-a-Year Control","Location","nw")

%%
function output = dydt_uncontrolled(t,y,vars)
    x = y(1);
    dx = y(2);
    x_target = vars(8);
    if t > 5 
        kp = 0;
        kd = 0;
        ki = 0;
    else
        kp = 0;
        kd = 0;
        ki = 0;
    end
    sinusodal_forcing = -vars(7)*(vars(2))^2*sin(vars(2)*t);
    f = vars(1);
    d2x = -kp*(x-x_target)+sinusodal_forcing+f;
    output = [dx;d2x];
end

function output = dydt(t,y,vars)
    x = y(1);
    dx = y(2);
    x_target = vars(8);
    if t > 5 
        kp = vars(3);
        kd = vars(5);
    else
        kp = 0;
        kd = 0;
    end
    sinusodal_forcing = -vars(7)*(vars(2))^2*sin(vars(2)*t);
    f = vars(1);
    if t == 0
        t = realmin;
        
    end
    control_forcing = -kp*(x-x_target)-kd*(dx);
    d2x = control_forcing+sinusodal_forcing+f;
    output = [dx;d2x];
end

function output = dydt_sparse(t,y,vars)
    x = y(1);
    dx = y(2);
    x_target = vars(8);
    if t > 5 
        kp = vars(3);
        kd = vars(5);
    else
        kp = 0;
        kd = 0;
    end
    vars(6);
    sinusodal_forcing = -vars(7)*(vars(2))^2*sin(vars(2)*t);
    f = vars(1);
    control_forcing = -kp*(vars(4)-x_target)-kd*vars(6);
    d2x = control_forcing+sinusodal_forcing+f;
    output = [dx;d2x];

end