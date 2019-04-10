%% houskeeping

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

Mwater_i = (VWaterInit*RhoWater);
%% velocity of the wind as initial conditions

% x y z respectively:
%x = downrange.
%z = height.
%y = cross range.
Vwx = 0 ;
Vwy = 10 ;
Vwz = 0 ;

%initial conditions for ode:
VelX0 = 0;
VelZ0 = 0;
VelY0 = 0;

x0 = 0;
z0 = z0;
y0 = 0; %intial condition for location into the page


%% read data:

data = load('StaticTestData/LA8am_test6') .* 4.44822   ; % to convert to newton
time_test = linspace(0,length(data)/1652,length(data));

% delte the noise and not useful data:

figure(1)

plot(time_test,data(:,3),'r.-','LineWidth',0.5)
[ x y ] = ginput(2);
[ ix1 iy1 ] = min(abs(time_test - x(1)));
[ ix2 iy2 ] = min(abs(time_test - x(2)));



Thrust_test = data(iy1:iy2,3);
time_test = time_test(iy1:iy2);
% zero time
time_test = time_test - time_test(1) ;

im = 0.00012;

% %find first min:
% difference = abs(time_test-im);
% ind1 = find(difference==min(difference));
% % zero and find the second min:
% difference(ind1); = [];
% ind2 = find(difference==min(difference));
% 
% Thrust = (Thrust_test(ind1) + Thrust_test(ind) )/2;

Opts = odeset('Events',@HitGround);
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThrustODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz,Thrust_test,time_test), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);


%% plot height and range:

figure(2);
plot(Results(:,7),Results(:,8),'-','Color',[0.25 0.25 0.25],'LineWidth',1.4)
hold on

MaxHeight = find(Results(:,8)==max(Results(:,8)),7);
plot(Results(MaxHeight,7),max(Results(:,8)),'o','Color',[ 0 0.5 0],'MarkerSize',7,'MarkerFaceColor',[0 0.5 0])
text(Results(MaxHeight,7)+0.5,max(Results(:,8))+0.5, ['Max Height =' num2str(max(Results(:,8))) 'm']);

plot(Results(end,7),0,'o','Color',[1 0 0],'MarkerSize',7,'MarkerFaceColor',[1 0 0])
text(Results(end,7)-21,0.5, ['Max Range =' num2str(Results(end,7)) 'm']);
ylim([0 floor(max(Results(:,8)))+3])
xlim([0 floor(Results(end,7))+3])
grid minor
title('Range vs Height')
xlabel('Range (m)')
ylabel('Height (m)')
legend('Rocket trajectory','Max Height','Max Range','Location','NorthWest')

%% consider 3d plotting after adjusting velocities of wind and made it 3d

figure(3) ; 

plot3(Results(:,7),Results(:,9),Results(:,8),'-','Color',[0.25 0.25 0.25],'LineWidth',3)
hold on
cord1 = [ -5 -50 0 ];
cord2 = [ -5 50 0 ];
cord3 = [ 100 50 0 ];
cord4 = [ 100 -50 0 ];

img = imread('CUField.jpg');

points = [ cord1' cord2' cord3' cord4' ];
%fill3(points(1,:),points(2,:),points(3,:),'g')
alpha(0.2)
hold on
%image('CData',img,'XData',[cord1(1) cord2(1) cord3(1) cord4(1)],'YData',[cord1(3) cord2(3) cord3(3) cord4(3)])

surf([0 109.1; 0 109.1],[24.25 24.25; -24.25 -24.25]...
    ,[0 0; 0 0],...    % Plot the surface
     'CData',img,...
     'FaceColor','texturemap');

grid on


xlabel('Downrange(m)')
ylabel('Crossrange(m)')
zlabel('Height (m)')
view([-30 45]) 

title('3D Flight path')

