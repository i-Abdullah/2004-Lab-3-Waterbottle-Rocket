%% info:


%{

this code will run uncertainty analysis using monter carlo simulations.

%}

%% housekeeping

clear
clc
close all;

%% ODE conditions: Launch 1


%% Launch 1:

% {
g = 9.81; % m/s2, acceleration due to gravity,
Cd= 0.8; % discharge coefficient
Rhoairamb = 0.961; % kg/m^3 ambient air density
Volbottle= 0.002; % m^3 volume of empty bottle
Pamb= 12.1*6894.76; % converted to Pa, atmospheric pressure
GammaGas = 1.4; % ratio of specific heats for air
RhoWater = 1000; % kg/m^3, density of water
DThroat= 2.1; % cm, diameter of throat
DBottle= 10.5; % in cm, diameter of bottle
R = 287; %J/kgK, gas constant of air
MBottle= 0.15; % kg mass of empty 2-liter bottle with cone and fins
MBottle= 0.110; % kg mass of empty 2-liter bottle with cone and fins
CD= 0.35; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.0006; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
TAirInit = 290.15; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 45 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
z0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
PayLoad = 0 ;
Fins = 0 ;
TotalMass0 = PayLoad + Fins + MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

TotalMass0 = 694e-3 ;
%MassAirInit = TotalMass0 - 600e-3 - MBottle  ;

% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

% the following from TA's Launch

WindSpeed = 9; %mph
[ xwind ywind ] = WindLaunch(204, 'ENE',WindSpeed);
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;

%initial conditions for ode:

VelX0 = 0;
VelZ0 = 0;
VelY0 = 0;

x0 = 0;
z0 = z0;
y0 = 0; %intial condition for location into the page


%}


%% get actual landing prediction:

Opts = odeset('Events',@HitGround);
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Estimatedx = Results(end,7);
Estimatedy = Results(end,9);

% actual landing from the day of the flight

Actualx = 65;
Actualy = Actualx * tand(7) ; % 7 is the drift angle, measured the same day.


Correctedx = Actualx;
Correctedy = Estimatedy ; % 7 is the drift angle, measured the same day.

%% redefine the terms that has uncertintiy

N = 1000;

CD = randn(N,1)*0.05 + CD ;
Pgage = randn(N,1)*2*6894.76 + Pgage ;
VWaterInit = randn(N,1)*0.00001 + VWaterInit ;
Theta = randn(N,1)*1 + Theta ;
WindSpeed = randn(N,1)*1.5 + WindSpeed;

% Now randoamly sample even the 

PossibleWind = { 'ENE' ; 'NE' ; 'NNE' };

%% run simulations
for i=1:N
    
VAirInit = Volbottle - VWaterInit(i);
 
[ xwind ywind ] = WindLaunch(204, char(PossibleWind(randi(numel(PossibleWind)))),WindSpeed(i));
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;



Opts = odeset('Events',@HitGround);
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta(i),Pgage(i),Pamb,Cd,ThroatArea,CD(i),BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

% Store results:

x(i) = Results(end,7);
y(i) = Results(end,9);

    
end


%% make error ellipses:

subplot(2,1,1)

plot(x,y,'k.','markersize',5)
axis equal;
grid on;
xlabel('x (Downrange) (m)');
ylabel('y (Crossrange) (m)');

hold on;

% Calculate covariance matrix

P = cov(x,y);
mean_x= mean(x);
mean_y= mean(y);% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
[eigvec,eigval] = eig(P); % Compute eigen-stuff

xy_vect= [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect= xy_vect(:,1);y_vect= xy_vect(:,2);% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b-.','LineWidth',2)
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g-.','LineWidth',2)
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r-.','LineWidth',2)

plot(Estimatedx,Estimatedy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[1, 1, 0],'MarkerSize',9)

plot(Actualx,Actualy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[0.3 0.8 0],'MarkerSize',9)

plot(Correctedx,Correctedy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[1 0.2 0],'MarkerSize',9)

grid minor

legend('Simulated data','1 \sigma','2 \sigma','3 \sigma','Estimated landing position','Actual landing position',...
    'Corrected landing position','Location','SouthEast')
title('Error ellipses for landing position') 

subplot(2,1,2)

plot(x,y,'k.','markersize',5)
axis equal;
grid on;
xlabel('x (Downrange) (m)');
ylabel('y (Crossrange) (m)');

hold on;

% Calculate covariance matrix

P = cov(x,y);
mean_x= mean(x);
mean_y= mean(y);% Calculate the define the error ellipses
n=100; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
[eigvec,eigval] = eig(P); % Compute eigen-stuff

xy_vect= [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect= xy_vect(:,1);y_vect= xy_vect(:,2);% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b-.','LineWidth',2)
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g-.','LineWidth',2)
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r-.','LineWidth',2)

plot(Estimatedx,Estimatedy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[1, 1, 0],'MarkerSize',9)

plot(Actualx,Actualy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[0.3 0.8 0],'MarkerSize',9)

plot(Correctedx,Correctedy,'mo','MarkerEdgeColor','k',...
'MarkerFaceColor',[1 0.2 0],'MarkerSize',9)

grid minor

title('Zoomed-in plot') 
ylim([ 0.5 2.5])
