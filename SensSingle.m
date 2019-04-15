%% info:

% Date Created : 24 Nov, 2018.

%{
- Done by:

%}

% this scripts is just for the purpose of plotting the parameters of water bottle
% rockets and how it affects max height and range. for more info read
% project2.m

clear;
clc;
close all;


%% inital conditions 


% updated conditions:

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
CD= 0.5; % drag coefficient
Pgage= 50*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
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

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

Vwx = 0 ;
Vwy = 0 ;
Vwz = 0 ;



%% 

%initial conditions for ode:

VelX0 = 0;
VelZ0 = 0;
VelY0 = 0;

x0 = 0;
z0 = z0;
y0 = 0; %intial condition for location into the page

Opts = odeset('Events',@HitGround);


%% How angel of launch affects max height and range

%run ODE45 for all values of theta 1 to 90
for i=1:90
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,i,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Var(i) = i;
Xmax(i) = Results(end,7); % Max Range;
Zmax(i) = max(Results(:,8)); % maximum height

end


% find maximum X_location:

%plot max height and max distance on same plot for each theta
figure(1);

plot(Var,Xmax,'LineWidth',2)
hold on;
plot(Var,Zmax,'LineWidth',2)
hold on
plot(Theta,Xmax(find(Var==Theta)),'or','MarkerSize',5,'MarkerFaceColor','r')
plot(Var(find(max(Xmax)==Xmax)),Xmax(find(max(Xmax)==Xmax)),'go','MarkerSize',5,'MarkerFaceColor','g')
line([10 10],[0 70],'LineWidth',2,'Color','k','LineStyle',':');
line([70 70],[0 70],'LineWidth',2,'Color','k','LineStyle',':');

legend('Max Range', 'Max Height',['Baseline=' num2str(Theta) char(176)],['Optimal=' num2str(Var(find(max(Xmax)==Xmax))) char(176)],...
    'Lowe bound','Upper bound');
title('\theta vs Max Height and Range');
xlabel('\theta (degrees)');
ylabel('Distance (meters)');
grid minor


% store derivatives to see how the downrange and height changes with
% respect to the considered variable

dtheta = Var; % store the variables
dx_dtheta = diff(Xmax); % store X derivatives.
dz_dtheta = diff(Zmax); % Store Z derivatives


%% How coefficient of drag affects max height and range

% clear variable:
clear Zmax
clear Xmax
clear Var


%ode initial condition

%run ODE45 for all values of CD = 0.1 to 1;
i = 1; % index
for j=0.1:0.01:1
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,j,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Var(i) = j;
Xmax(i) = Results(end,7); % Max Range;
Zmax(i) = max(Results(:,8)); % maximum height
i = i+1;
end


%plot max height and max distance on same plot for each theta
figure(2)

plot(Var,Xmax,'LineWidth',2)
hold on;
plot(Var,Zmax,'LineWidth',2)
plot(CD,Xmax(find(Var==CD)),'or','MarkerSize',5,'MarkerFaceColor','r')
plot(Var(find(max(Xmax)==Xmax)),Xmax(find(max(Xmax)==Xmax)),'go','MarkerSize',5,'MarkerFaceColor','g')

line([0.2 0.2],[0 120],'LineWidth',2,'Color','k','LineStyle',':');
line([0.9 0.9],[0 120],'LineWidth',2,'Color','k','LineStyle',':');


legend('Max Range', 'Max Height',['Baseline=' num2str(CD)],['Optimal=' num2str(Var(find(max(Xmax)==Xmax)))],...
    'Lowe bound','Upper bound');
title('Coefficient of Drag vs Max Height and Range');
xlabel('Coefficient of Drag');
ylabel('Distance (meters)');
grid minor




%% How volume of water affects max height and range


% clear variable:
clear Zmax
clear Xmax
clear Var


%ode initial condition

%run ODE45 for all values of Volume of water of 0.0001 m^3 to 0.002;
i = 1; % index

% WHEN VOLUME INCREASED TO MAX VOLUME OF BOTTLE IT IS NOT WORKING!!
for j=.0001:.00001:.0013
% Call ODE
    VAirInit=Volbottle-j;
    TotalMass0 = MBottle + (j*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)) ; % Total mass
    MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));
j
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 40],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Var(i) = j;
Xmax(i) = Results(end,7); % Max Range;
Zmax(i) = max(Results(:,8)); % maximum height
i = i+1;
end


%reset the changing variable

VWaterInit= 0.001; % m^3, initial volume of water inside bottle
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air


%plot max height and max distance on same plot for each theta
figure(3)
plot(Var,Xmax,'LineWidth',2)
hold on;
plot(Var,Zmax,'LineWidth',2)

plot(VWaterInit,Xmax(find(Var==VWaterInit)),'or','MarkerSize',5,'MarkerFaceColor','r')
plot(Var(find(max(Xmax)==Xmax)),Xmax(find(max(Xmax)==Xmax)),'go','MarkerSize',5,'MarkerFaceColor','g')

line([1e-4 1e-4],[0 120],'LineWidth',2,'Color','k','LineStyle',':');
line([1e-3 1e-3],[0 120],'LineWidth',2,'Color','k','LineStyle',':');


legend('Max Range', 'Max Height',['Baseline=' num2str(VWaterInit)],['Optimal=' num2str(Var(find(max(Xmax)==Xmax)))],...
    'Lowe bound','Upper bound');

title('Volume of water vs Max Height and Range');
xlabel('Volume of water (m^3)');
ylabel('Distance (meters)');
grid minor





%% How gage pressure affects max height and range


% clear variable:
clear Zmax
clear Xmax
clear Var


%ode initial condition

%run ODE45 for all values of Volume of water of 0.0001 m^3 to 0.002;
i = 1; % index
for j=25:1:70
% Call ODE
    Pgage= j*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
    TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
    MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));

[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Var(i) = Pgage;
Xmax(i) = Results(end,7); % Max Range;
Zmax(i) = max(Results(:,8)); % maximum height
i = i+1;
end

%Reset variables
Pgage=50*6894.76;
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air



%plot max height and max distance on same plot for each theta
figure(4)

plot(Var,Xmax,'LineWidth',2)
hold on;
plot(Var,Zmax,'LineWidth',2)

plot(Pgage,Xmax(find(Var==Pgage)),'or','MarkerSize',5,'MarkerFaceColor','r')
plot(Var(find(max(Xmax)==Xmax)),Xmax(find(max(Xmax)==Xmax)),'go','MarkerSize',5,'MarkerFaceColor','g')

line([2e5 2e5],[0 80],'LineWidth',2,'Color','k','LineStyle',':');
line([4e5 4e5],[0 80],'LineWidth',2,'Color','k','LineStyle',':');


legend('Max Range', 'Max Height',['Baseline=' num2str(Pgage)],['Optimal=' num2str(Var(find(max(Xmax)==Xmax)))],...
    'Lowe bound','Upper bound');

title('Gage pressure vs Max Height and Range');

xlabel('Gage pressure (Pa)');
ylabel('Distance (meters)');

grid minor


%% How change in temp of fuel changes the? 

% clear variable:
clear Zmax
clear Xmax
clear Var


%ode initial condition

%run ODE45 for all values of Volume of water of 0.0001 m^3 to 0.002;
i = 1; % index
for j=283:1:313
% Call ODE
TAirInit = j;


MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));


[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);

Var(i) = TAirInit;
Xmax(i) = Results(end,7); % Max Range;
Zmax(i) = max(Results(:,8)); % maximum height
i = i+1;
end

%Reset variables
TAirInit = 300;
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air



%plot max height and max distance on same plot for each theta
figure(5)

plot(Var,Xmax,'LineWidth',2)
hold on;
plot(Var,Zmax,'LineWidth',2)

plot(TAirInit,Xmax(find(Var==TAirInit)),'or','MarkerSize',5,'MarkerFaceColor','r')
plot(Var(find(max(Xmax)==Xmax)),Xmax(find(max(Xmax)==Xmax)),'go','MarkerSize',5,'MarkerFaceColor','g')

line([284.5 284.5],[0 80],'LineWidth',2,'Color','k','LineStyle',':');
line([310 310],[0 80],'LineWidth',2,'Color','k','LineStyle',':');


legend('Max Range', 'Max Height',['Baseline=' num2str(TAirInit)],['Optimal=' num2str(Var(find(max(Xmax)==Xmax)))],...
    'Lowe bound','Upper bound');

title('Water Temp vs Max Height and Range');

xlabel('Water Temp (K)');
ylabel('Distance (meters)');

grid minor
