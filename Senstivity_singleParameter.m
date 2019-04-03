%% info:

% Date Created : 24 Nov, 2018.

%{
- Done by:
            1- Brendan Palmer, id : 108102169
            2- Abdulla AlAmeri id : 109364560

%}

% this scripts is just for the purpose of plotting the parameters of water bottle
% rockets and how it affects max height and range. for more info read
% project2.m

clear;
clc;
close all;


%% inital conditions 


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
Pgage= 66.5*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 43.1 ; % initial angle of rocket in degress
X0 = 0.0; % in meters, initial horizontal distance
z0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
PayLoad = 25*10^-3 ;
Fins = 10*10^-3 ;
TotalMass0 = PayLoad + Fins + MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air


%% Numerical integration.

%% How angel of launch affects max height and range

%ode initial condition

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = z0;
x=[];
ymax=[];
xmax=[];

%run ODE45 for all values of theta 1 to 90
for i=1:90
    [ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,i,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 6],[TotalMass0 MassAirInit...
    VAirInit VelX0 VelZ0 Range0 z0 ]);
    %add the values for theta into the x array
    x=cat(1,x,i);
    %find max height and add to y array
    ymax=cat(1,ymax,max(Results(:,7)));
    %initialize j
    j=1;
    %find where y value is first negative
    while Results(j,7)>0
        j=j+1;
    end
    %interpolate max distance
    m=(Results(j,7)-Results(j-1,7))/(Results(j,6)-Results(j-1,6));
    %concatinate the max distance with max distance for specific theta
    xmax=cat(1,xmax,Results(j-1,6)+(0-Results(j-1,7))/m);
end
%plot max height and max distance on same plot for each theta
figure;
subplot(2,2,1)
plot(x,ymax)
hold on;
plot(x,xmax)
legend('Max height', 'Max Range');
title('\theta vs Max Height and Range');
xlabel('\theta (degrees)');
ylabel('Distance (meters)');
grid minor
%reset the changing variable
Theta=45;

%% How coefficient of drag affects max height and range

%ode initial condition

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;
x=[];
ymax=[];
xmax=[];
%run ODE45 for all values of drag coefficient to 0 to 1
for i=0:.001:1
    [ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,i,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 6],[TotalMass0 MassAirInit...
    VAirInit VelX0 VelZ0 Range0 y0 ]);
    %add the values for drag coefficient into the x array
    x=cat(1,x,i);
    %find max height and add to y array
    ymax=cat(1,ymax,max(Results(:,7)));
    %initialize j
    j=1;
    %find where y value is first negative
    while Results(j,7)>0
        j=j+1;
    end
    %interpolate max distance
    m=(Results(j,7)-Results(j-1,7))/(Results(j,6)-Results(j-1,6));
    %concatinate the max distance with max distance for specific CD
    xmax=cat(1,xmax,Results(j-1,6)+(0-Results(j-1,7))/m);
end
%plot max height and max distance on same plot for each CD
%figure;
subplot(2,2,2)
plot(x,ymax)
hold on;
plot(x,xmax)
legend('Max height', 'Max Range');
title('Cefficient of Drag vs Max Height and Range');
xlabel('Drag Coefficient');
ylabel('Distance (meters)');
grid minor
CD=.5;

%% How volume of water affects max height and range

%ode initial condition

VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;
x=[];
ymax=[];
xmax=[];
%run ODE45 for all values of volume water to 0 to .002
for i=.0001:.001:.002
    %use i to find the volume of air in the bottle
    VAirInit=Volbottle-i;
    TotalMass0 = MBottle + (i*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
    MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit));
    [ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 6],[TotalMass0 MassAirInit VAirInit VelX0 VelZ0 Range0 y0 ]);
    %add the values for drag coefficient into the x array
    x=cat(1,x,i);
    %find max height and add to y array
    ymax=cat(1,ymax,max(Results(:,7)));
    %initialize j
    j=1;
    %find where y value is first negative
    while Results(j,7)>0
        j=j+1;
    end
    %interpolate max distance
    m=(Results(j,7)-Results(j-1,7))/(Results(j,6)-Results(j-1,6));
    %concatinate the max distance with max distance for specific CD
    xmax=cat(1,xmax,Results(j-1,6)+(0-Results(j-1,7))/m);
end
%plot max height and max distance on same plot for each CD
%figure;
subplot(2,2,3)
plot(x,ymax)
hold on;
plot(x,xmax)
legend('Max Height', 'Max Range');
title('Initial Volume Water vs Max Height and Range');
xlabel('Volume Water Initial (m^3)');
ylabel('Distance (meters)');
grid minor
VWaterInit= 0.001; % m^3, initial volume of water inside bottle
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

%% How gage pressure affects max height and range

%ode initial condition

y0=.25;
VelX0 = 0;
VelZ0 = 0;
Range0 = 0;
Height0 = y0;
x=[];
ymax=[];
xmax=[];
%run ODE45 for all values of gage pressure from 25 to 75 psi
for i=25:75
    Pgage= i*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
    TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((i+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
    MassAirInit = (((i+Pamb)*VAirInit ) / (R*TAirInit));
    [ Time Results ] = ode45(@(Time,States) RocketODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,y0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R), [ 0 5 ],[TotalMass0 MassAirInit...
    VAirInit VelX0 VelZ0 Range0 y0 ]);
    %add the values for drag coefficient into the x array
    x=cat(1,x,i);
    %find max height and add to y array
    ymax=cat(1,ymax,max(Results(:,7)));
    %initialize j
    j=1;
    %find where y value is first negative
    while Results(j,7)>0
        j=j+1;
    end
    %interpolate max distance
    m=(Results(j,7)-Results(j-1,7))/(Results(j,6)-Results(j-1,6));
    %concatinate the max distance with max distance for specific CD
    xmax=cat(1,xmax,Results(j-1,6)+(0-Results(j-1,7))/m);
end
%plot max height and max distance on same plot for each CD
%figure;
subplot(2,2,4)
plot(x,ymax)
hold on;
plot(x,xmax)
legend('Max Height', 'Max Range');
title('Pressure vs Max Height and Range');
xlabel('Pressure (Pa)');
ylabel('Distance (meters)');
grid minor
%Reset variables
Pgage=50*6894.76;
TotalMass0 = MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air
