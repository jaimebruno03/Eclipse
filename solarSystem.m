function [r,body] = solarSystem( t1, t2, dt )

% Generates trajectories for selected bodies of the solar system.
% J.M.Soler. Dec.2016
% 
% Input:
% t1  : initial time of trajectories, from 0h, Jan.1, 2017, in days
% t2  : final time of trajectories, from 0h, Jan.1, 2017, in days
% dt  : integration time interval, in days
%
% Output:
% r(3,nb,nt) : vector positions of nb bodies at nt times, in au
% body(nb,9) : name of the nb bodies, with one trailing blank
%
% Algorithm:
% Positions 3 hours before and after t0 = 3h, Jan.1, 2017, are taken from
%   http://aa.usno.navy.mil/data/docs/geocentric.php  (Unavailable since 2019)
%   https://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/sun_en.cgi    (Sun)
%   https://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/planet_en.cgi (Planets)
%   https://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/moon_en.cgi   (Moon)
% Their mean and difference are used to compute r0 and v0.
% Then we integrate Newton's equations, using Verlet method, 
% from t0 to t1 and we find r1 and v1 at t1. Finally, we find the 
% trajectory from t1 to t2 with the same method.

% Set auxiliary variables and arrays
day = 24*60*60;             % one day in s
km = 1.e3;                  % one km in m
au = 1.495978707e11;        % one astronomical unit, in m
eclipticObliquity = 23.44;  % angle between Earth orbit and equator, in deg

% Set initial positions at three times, relative to the Earth center.
% mass: in kg
% ascension: angle along the equator, from the equinox, in h,min,sec
% declination: angle from the equator, in deg,min,sec
% distance: in astronomical units (au)
t0 = 3*60*60;    % time of data, in s (3h, Jan.1, 2017)
dt0 = 3*60*60;   % 3h time interval between the data positions, in s

nb = 1; body(nb,:) = sprintf('%8s ','Sun');
mass(nb) = 1.9891e30;
ascension(:,:,nb) = [18 46 46.974      % position at t0-dt0
                     18 47 20.092      % position at t0
                     18 47 53.205];    % position at t0+dt0
declination(:,:,nb) = -[22 59 56.65
                        22 59 19.47
                        22 58 41.85];
distance(:,nb) = [0.983337917
                  0.983336025
                  0.983334194]*au;

nb = nb+1; body(nb,:) = sprintf('%8s ','Earth');
mass(nb) = 5.9736e24;
ascension(:,:,nb) = [0 0 0     % position zero, since Earth is the origin
                     0 0 0
                     0 0 0];
declination(:,:,nb) = [0 0 0
                       0 0 0
                       0 0 0];
distance(:,nb) = [0
                  0
                  0];

%nb = nb+1; body(nb,:) = sprintf('%8s ','Mars');
%mass(nb) = 6.4185e23; 
%ascension(:,:,nb) = [22 45 41.412
%                     22 46 02.449
%                     22 46 23.483];
%declination(:,:,nb) = -[8 48 35.59
%                        8 46 20.29
%                        8 44 04.93];
%distance(:,nb) = [1.640451805
%                  1.641286008
%                  1.642120297]*au;
%              
% nb = nb+1; body(nb,:) = sprintf('%8s ','Venus');
% mass(nb) = 4.868e24; 
% ascension(:,:,nb) = [
%                      
%                      ];
% declination(:,:,nb) = -[
%                         
%                         ];
% distance(:,nb) = [
%                   
%                   ]*au;
% 
 nb = nb+1; body(nb,:) = sprintf('%8s ','Moon');
 mass(nb) = 7.3477e22; 
 ascension(:,:,nb) = [20 54 53.859
                      21 01 22.256
                      21 07 50.227];
 declination(:,:,nb) = -[15 20 18.36
                         15 01 17.64
                         14 41 33.72];
 distance(:,nb) = [391301
                   390880
                   390457]*km;

% Find initial positions and velocities as 3D vectors
for ib = 1:nb
    phi = ascension(:,:,ib) * [60*60; 60; 1] * (2*pi)/(24*60*60);  % in rad
    theta = declination(:,:,ib) * [60*60; 60; 1] * (2*pi)/(360*60*60) ...
          + pi/2;      % angle from the north pole, in rad
    rd = zeros(3,3);   % relative coordinates at three times
    rd(:,1) = distance(:,ib) .* sin(theta) .* cos(phi);  % x coordinates
    rd(:,2) = distance(:,ib) .* sin(theta) .* sin(phi);  % y
    rd(:,3) = distance(:,ib) .* cos(theta);              % z
    rd = rd';                              % change row to column vectors
    r0(:,ib) = rd(:,2);                    % position vectors at t0
    v0(:,ib) = (rd(:,3)-rd(:,1))/(2*dt0);  % velocity vectors at t0
end

% Rotate around equinox (x axis) to place Earth's orbit in xy plane
% angle = eclipticObliquity * pi/180;
% rotMat = [1     0          0
%           0 cos(angle) -sin(angle)
%           0 sin(angle)  cos(angle)];
% r0 = rotMat*r0;
% v0 = rotMat*v0;

% Change to center of mass system
rcm = r0*mass' / sum(mass);           % center of mass position
vcm = v0*mass' / sum(mass);           % center of mass velocity
for ib = 1:nb
    r0(:,ib) = r0(:,ib) - rcm;
    v0(:,ib) = v0(:,ib) - vcm;
end

% Change input times from days to seconds
t1 = t1*day;
t2 = t2*day;
dt = dt*day;

% Find position and velocity at t1
myForce = @(r)gravity(mass,r);
if (t1==t0)
    
    r1 = r0;
    v1 = v0;

else
    
    % Find trajectories from t0 to t1
    nt = max(1,round(abs(t1-t0)/dt));  % number of time intervals up to t1
    dt1 = (t1-t0)/nt;                  % exact time interval
    [r,v] = verlet( myForce, mass, r0, v0, dt1, t1-t0 );
    
    % Find position and velocity at t1
    r1 = r(:,:,end);             % position at t1
    v1 = v(:,:,end);             % velocity at t1

end

% Find trajectory between t1 and t2
r = verlet( myForce, mass, r1, v1, dt, t2-t1 );
r = r / au;                          % Change from m to au

end % function solarSystem
