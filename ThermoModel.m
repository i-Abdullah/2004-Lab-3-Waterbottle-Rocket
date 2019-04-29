    %% info

%{
% This's the main script that will be using numerical inegration to model
water bottle rocket.

% This code is using ode45 to analyze the motion of water bottle rocket, it
uses initial conditions and run numerical integration via ode45 to analyze the
rocket.

Purpose: Discover how to use derivatives and numerical integration to
determine the rocket thrust and trajectory. Learn how rocket specifications
and change trajectory.

Givens: Needed x displacement of the rocket. Derivatives modeling change in
mass, pressure, volume, and velocity

Required: Determine the launch angle, and bottle rocket dimensions to
reach a particular distance.

Assumptions: Expansion of air is isentropic, no change of rocket deflection
(in y- direction), rocket adjusts flight path to match the heading
instantaneously, and static air around the rocket.


ASEN 2012, Fall 18.

- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

For ASEN 2004: Vehicle Design:

3D simulations are included. 


%}

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

%% Launch 1:

%{
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

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

% the following from TA's Launch

[ xwind ywind ] = WindLaunch(204, 'ENE',9);
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;

%}

%% Launch2 :

%{
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
CD= 0.324; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.000594 ; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
TAirInit = 290.15; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 43 ; % initial angle of rocket in degress
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

TotalMass0 = 704e-3 ;
%MassAirInit = TotalMass0 - 600e-3 - MBottle  ;

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

% the following from TA's Launch

[ xwind ywind ] = WindLaunch(204, 'SSW',2);
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;

%}


%% Launch3 :

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
CD= 0.324; % drag coefficient
Pgage= 40*6894.76; % in pascal, the 6894.76 is to convert. initial gage pressure of air in bottleVolwater,
VWaterInit= 0.001 ; % m^3, initial volume of water inside bottle
TAirInit = 300; % K, initial temperature of
TAirInit = 290.15; % K, initial temperature of
Airv0 = 0.0 ;% m/s, initial velocity of rocket
Theta= 35 ; % initial angle of rocket in degress
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

TotalMass0 = MBottle + 1 ;
%MassAirInit = TotalMass0 - 600e-3 - MBottle  ;

%% velocity of the wind as initial conditions

% x y z respectively:

%x = downrange.
%z = height.
%y = cross range.

% the following from TA's Launch

[ xwind ywind ] = WindLaunch(204, 'SW',3);
Vwx = xwind*0.44704 ; %multiplication to convert from mph to m/s
Vwy = ywind*0.44704 ;
Vwz = 0 ;

%}
%% 

%initial conditions for ode:

VelX0 = 0;
VelZ0 = 0;
VelY0 = 0;

x0 = 0;
z0 = z0;
y0 = 0; %intial condition for location into the page

Opts = odeset('Events',@HitGround);
% Call ODE
[ Time Results ] = ode45(@(Time,States) ThermoODE(Time,States,TestStandLength,Theta,Pgage,Pamb,Cd,ThroatArea,CD,BottleArea,Rhoairamb,RhoWater,Volbottle,z0,VAirInit,GammaGas,g,TAirInit,MassAirInit,R,Vwx,Vwy,Vwz), [ 0 5],[TotalMass0 MassAirInit...
VAirInit VelX0 VelZ0 VelY0 x0 z0 y0 ],Opts);



%% calculate thrust

%ODE Has really a hard way of calculating the time, even though matrices t1
%and t2 can tell us the end time of each phase there's no time that matches
%it exactly in the Time matrix, and extracting the force of thrust via
%globals really yields bad results, thus the forces has to be calculated
%manually.


%index namings: those will use to index the storing of thrust and time
store1 = 0;
store2 = 0;
store3 = 0;

[ r c ] = size(Results);
for i=1:r
    
if Results(i,3) < Volbottle
    
Pressure1(i) = ( ( VAirInit ./ Results(i,3) ) .^ GammaGas ) .* (Pgage+Pamb) ; 
Thrust1(i) = 2.* Cd .* ThroatArea .* ( Pressure1(i) - Pamb) ;
TP1(i) = Time(i);

%phase 2
elseif Results(i,3)>= Volbottle
    
    %T and P of end states
Tend = TAirInit * (( VAirInit/Volbottle) ^ (GammaGas-1) );
Pend = (Pgage+Pamb) * (( VAirInit/Volbottle) ^ (GammaGas) );
    
PressureCond = Pend * (Results(i,2)/MassAirInit)^(GammaGas) ;

if PressureCond>Pamb
Density = Results(i,2) / Volbottle;
Temp = PressureCond/(Density*R);
CriticalP = (PressureCond) * (2./(GammaGas+1)).^(GammaGas/(GammaGas-1));
  
if CriticalP > Pamb
    
    Mach  = 1;
    Texit = (2/(GammaGas+1))*Temp ;
    Vexit = sqrt(GammaGas*Texit*R);
    Pexit = CriticalP;
    Densityexit = CriticalP/(R*Texit) ;
    
elseif CriticalP <= Pamb
    
   Mach = sqrt(( (PressureCond/Pamb)^( ( (GammaGas-1)/GammaGas)) - 1 ) * (2/(GammaGas-1)));
   Texit = Temp/(1+((GammaGas-1)/2)*Mach^2);
   Pexit = Pamb;
   Densityexit = Pamb/(R*Texit) ;
   Vexit = Mach * sqrt(GammaGas*Texit*R);
end

% how mass of air and rocket change with time

MassAirFlowRate = Cd*Densityexit*ThroatArea*Vexit;

store2 = store2 + 1;
Thrust2(store2) = MassAirFlowRate *Vexit + (Pexit-Pamb)*ThroatArea ;
TP2(store2) = Time(i);

else
store3 = store3 + 1;
Thrust3(store3) = 0 ;
TP3(store3) = Time(i);


%% Phase 3: 
end




end

end


%% plot thrust:

figure(1);
plot([TP1 TP2 TP3],[Thrust1 Thrust2 Thrust3],'Color',[0.25 0.25 0.25],'LineWidth',1.4)
hold on
plot(TP1(end),Thrust1(end),'o','Color',[ 0 0.5 0],'MarkerSize',7,'MarkerFaceColor',[0 0.5 0])
plot(TP2(end),Thrust2(end),'o','Color',[1 0 0],'MarkerSize',7,'MarkerFaceColor',[1 0 0])
xlim([0 TP3(floor(length(TP3)/4))])
grid minor
title('Thrust VS Time')
xlabel('Time (Seconds)')
ylabel('Thrust (N)')
legend('Thrust profile','End of Water Phase','End of Air Phase','Location','NorthEast')


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



%% for future devlopments, ignore.

%{
% this portion of the graph was taken from MIT matlab training section,

% http://web.mit.edu/8.13/matlab/MatlabTraining_IAP_2012/AGV/DemoFiles/ScriptFiles/html/Part3_Animation.html

% figure;
% 
% subplot(2,1,1);
% xlabel('time (sec)'); ylabel('angle (\circ)');
% 
% for id = 1:length(Time3)
%    subplot(2,1,1);
%    plot(Results(id,6),Results(id,7), 'LineWidth', 2);
%    line(Results(id,6), Results(id,7), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
%    line(Results2(id,6), Results2(id,7), 'Marker', '*', 'MarkerSize', 20, 'Color', 'r');
%    line(Results3(id,6), Results3(id,7), 'Marker', '^', 'MarkerSize', 20, 'Color', 'b');
% 
% %    line(Time(id), Results(id,7), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
%    xlabel('Range (m)'); ylabel('Height (m)');
% 
% %    % The bottom plot shows the animation of the double pendulum
% %    subplot(2,1,2);
% %    plot([0, x(id,1);x(id,1), x(id,2)], [0, y(id,1);y(id,1), y(id,2)], ...
% %       '.-', 'MarkerSize', 20, 'LineWidth', 2);
%    axis equal; axis([0 max(Results(:,6))+2 0 max(Results(:,7))+2]);
%    title(sprintf('Time: %0.2f sec', Time(id)));
%    grid minor
% 
%    drawnow;
%    
% end

%}
%{
clf;

ax = axes('XLim',[0 max(Results(:,6))+2],'YLim',[ 0 max(Results(:,7))+2],...
    'ZLim',[ 0 1 ]);

view(3)
grid minor;
axis equal

% create the body of the rocket:

% [ xcone ycone zcone ] = cylinder([0.2 -2]); % cylinder with top raduis of 0 so it's cone
% [ xcylin ycylin zcylin ] = cylinder([0.2 0.05]);

[ xcone ycone zcone ] = cylinder([0.4 0.05]); % cylinder with top raduis of 0 so it's cone
[ xcylin ycylin zcylin ] = cylinder([0.4 0.4]);

%create body parts of the rocket, assumption: no fins, we don't need them $_$

Parts(1) = surface(xcone,ycone,zcone,'FaceColor','red');
Parts(2) = surface(xcylin,ycylin,1.2*zcylin-1.2,'FaceColor','yellow');


%put the rocket together!
%SET BREAKPOINT HERE! ERROR IN THE NEXT LINE!!
Rocket = hgtransform('Parent',ax);
set(Parts,'Parent',Rocket)

%}