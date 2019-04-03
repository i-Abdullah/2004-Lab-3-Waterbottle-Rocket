%% info:

% Date Created : 20 March, 2019.

% this script attempts to model the trajectory of waterbottle rocket using 
% specific impulse (Isp) using static test stand data.
%% housekeeping

clear
clc
close all

%% define constants

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
y0 = 0.25; % in m, initial vertical height
TestStandLength= 0.5; % in m, length of test stand
VAirInit = Volbottle - VWaterInit ; %initial volume of Air.
ThroatArea = pi * ((DThroat*10^-2)/2)^2; %Area of throat
BottleArea  = pi * ((DBottle*10^-2)/2)^2; %Bottle Area
PayLoad = 0 ;
Fins = 0 ;
TotalMass0 = PayLoad + Fins + MBottle + (VWaterInit*RhoWater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air
Vwx = 1 ;
Vwy = 2 ;
Vwz = 1 ;


%% code:


data = load('StaticTestData/LA8am_test6') .* 4.44822   ; % to convert to newton

Mwater = VWaterInit*RhoWater;
%Mwater = 2.20462 ; 
%g = 32.2 ; 

%PosData = find(data(:,3)>0); %positive datap points

time = linspace(0,length(data)/1652,length(data));

sens_term = 100000; % senstivity term, adjust as needed.
[TF,S1,S2] = ischange(data(:,3),'linear','Threshold',sens_term);

i_end = find(TF==1); % index when there's big change at ends happend

% find where it started

sens_term = 10; % senstivity term, adjust as needed.
[TF,S1,S2] = ischange(data(:,3),'linear','Threshold',sens_term);

i_start = find(TF==1); % index when there's big change at ends happend


 time_indexed = time(i_start(1):i_end(end));
data_force = data(i_start(1):i_end(end),3);

ForceIntg = cumtrapz(time_indexed,data_force) ; %integral of force

Isp = ForceIntg(end) / (Mwater*g) ;
DV = Isp*g* log(TotalMass0/(PayLoad + Fins + MBottle)) ; % Change in velocity, assumes it's instantaneous.

%% ODE

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
MBottle= 117*10^-3; % kg mass of empty 2-liter bottle with cone and fins

CD= 0.5; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
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
Mwater = VWaterInit*RhoWater;
Mwater = 962*10^-3
Fins = 0 ;
TotalMass0 = PayLoad + Fins + MBottle + (Mwater) + (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); % Total mass
MassAirInit = (((Pgage+Pamb)*VAirInit ) / (R*TAirInit)); %initial mass of air

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

Vwx = 90 ;
Vwy = 0;
Vwz = 0 ;




%% 

%initial conditions for ode:

VelX0 = DV*cosd(Theta);
VelZ0 = DV*sind(Theta);
VelY0 = 0;

Range0 = 0;
Height0 = z0;
depth0 = 0; %intial condition for location into the page

%Stop = odeset('Events',HitGround);
% Call ODE
[ Time Results ] = ode45(@(Time,States) IspODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 Range0 z0 VelY0 depth0 ]);




%% make some plots :(

figure(1)

plot(time_indexed,data_force,'r.-','LineWidth',0.5)
hold on
area(time_indexed,data_force,'LineStyle',':','FaceAlpha',0.2,'FaceColor','r')

grid minor
title('Static test force of thrust')
ylabel('Time (seconds)')
xlabel('Thrust (N)')


figure(2);
plot(Results(:,6),Results(:,7),'-','Color',[1 0.5 0.2],'LineWidth',1.4)
hold on
LocationOFmax = Results((find(Results(:,7)==max(Results(:,7)))),6);
plot(LocationOFmax,max(Results(:,7)),'*','Color',[ 0 0.5 0],'MarkerSize',7,'MarkerFaceColor',[0 0.5 0])
text(LocationOFmax+0.5,max(Results(:,7))+0.5, ['Max Height =' num2str(max(Results(:,7))) 'm']);
j=1;
%find where y value is first negative
while Results(j,7)>0
    j=j+1;
end
%interpolate max distance
m=(Results(j,7)-Results(j-1,7))/(Results(j,6)-Results(j-1,6));
%concatinate the max distance with max distance for specific theta
xmax=Results(j-1,6)+(0-Results(j-1,7))/m;
plot(xmax,0,'*','Color',[0.2 0 0],'MarkerSize',7,'MarkerFaceColor',[0.2 0 0])
text(xmax-21,0.5, ['Max Range =' num2str(xmax) 'm']);
ylim([0 floor(max(Results(:,7)))+3])
xlim([0 floor(xmax)+3])
grid minor
title('Range vs Height')
xlabel('Range (m)')
ylabel('Height (m)')
legend('Rocket trajectory','Max Height','Max Range','Location','NorthWest')

MaxHeight = max(Results(:,7));
