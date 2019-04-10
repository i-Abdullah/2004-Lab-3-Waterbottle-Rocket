function [ x z ] = WindLaunch (theta, Direction,WindSpeed)
% this function will try to see how much of the wind is coming in the
% dirction of the launch pad for the waterbottle rocket, then convert it
% into components in a cartesian coordinate system with its origin at the
% launch pad.
%
% - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - 
% INPUTS:
%           - theta: angle of launch pad from Noth measured CCW to be +
%           - Direction: Direction of wind, is it N,S,W,E,NE,SW,..etc
%           - WindSpeed: Wind speed, output will be in the units of input
%
% - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - 
%
% Outputs:
%           - x z y : components of wind in a cartesian coordinate system
%                    with its origin at the launch pad.
%


% get unit vector along the launch pad:
% due to the weather station capabilities we can't get wind velocities
% in height direciton (z), we can only get downrange (x), crossrange (y).


hx = cosd(theta);
hy = sind(theta);



if theta <= 360
    
    switch Direction
        
        case 'N'
            x = hx*WindSpeed*cosd(0);
            z = hy*WindSpeed*sind(0);
            
        case 'E'
            x = hx*WindSpeed*cosd(90);
            z =hy*WindSpeed*sind(90);
            
        case 'S'
            x = hx*WindSpeed*cosd(180);
            z =hy*WindSpeed*sind(180);
            
        case 'W'
            
            x = hx*WindSpeed*cosd(270);
            z = hy*WindSpeed*sind(270);
            
        case 'NE'
            
            x = hx*WindSpeed*cosd(45);
            z = hy*WindSpeed*sind(45);
            
       case 'NNE'
            
            x = hx*WindSpeed*cosd(45/2);
            z = hy*WindSpeed*sind(45/2);
            
       case 'ENE'
            
            x = hx*WindSpeed*cosd(45 + 45/2);
            z = hy*WindSpeed*sind(45 + 45/2);
            
        case 'SE'
            
            x = hx*WindSpeed*cosd(90 + 45);
            z = hy*WindSpeed*sind(90 + 45);
            
        case 'ESE'
            
            x = hx*WindSpeed*cosd(90 + 45/2);
            z = hy*WindSpeed*sind(90 + 45/2);
            
        case 'SSE'
            
            x = hx*WindSpeed*cosd(90 + 45/2 + 45 );
            z = hy*WindSpeed*sind(90 + 45/2 + 45 );
            
            
        case 'SW'
            
            x = hx*WindSpeed*cosd(180 + 45 );
            z = hy*WindSpeed*sind(180 + 45 );
            
        case 'SSW'
            
            x = hx*WindSpeed*cosd(180 + 45/2 );
            z = hy*WindSpeed*sind(180 + 45/2 );

            
            
        case 'WSW'
            
            x = hx*WindSpeed*cosd(180 + 45 + 45/2 );
            z = hy*WindSpeed*sind(180 + 45 + 45/2 );
            
            
        case 'NW'
            
            x = hx*WindSpeed*cosd(270 + 45 );
            z = hy*WindSpeed*sind(270 + 45 );
            
            
        case 'NNW'
            
            x = hx*WindSpeed*cosd(270 + 45 + 45/2 );
            z = hy*WindSpeed*sind(270 + 45 + 45/2 );
            
            
        case 'WNW'
            
            
            x = hx*WindSpeed*cosd(270 + 45/2 );
            z = hy*WindSpeed*sind(270 + 45/2 );

     
    end
    
else
    
    errro('You can NOT put angle bigger than 360');
    
    
end

end


