function [value,isterminal,direction] = HitGround(Time,States)
% Locate the time when height passes through zero in a decreasing
% direction and stop integration

value = States(8); % detect when height (stored in 8th state) is =0;
isterminal= 1; % stop the integration
direction = -1; % negative direction