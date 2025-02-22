function [r,v,a] = verlet( force, mass, r0, v0, dt, tmax )

% Solves Newton's equation using velocity Verlet method for np particles 
% in nd dimensions. 
% 
% Input:
%   force     : function to calculate the force as f=force(r), in N
%   mass(1,np): particle masses of the np particles in kg
%   r0(nd,np) : intial positions of np particles in nd dimensions, in m
%   v0(nd,np) : initial velocities in m/s
%   dt        : integration time interval in s
%   tmax      : final integration time in s
% Output:
%   r(nd,np,nt) : position vectors at each time (0:dt:tmax), in m
%   v(nd,np,nt) : velocities at each time, in m/s
%   a(nd,np,nt) : accelerations at each time, in m/s^2
%                 If nd=1, the output sizes of x,v,a will be (np,nt)
%                 If np=1, the output sizes will be (nd,nt)
% J.M.Soler, Jan.2011

% Check array sizes
[nd,np] = size(r0);
if (numel(mass)~=np)
    error('verlet: ERROR: numel(mass)~=np')
end
if (size(v0)~=size(r0))
    error('verlet: ERROR: size(v0)~=size(x0)')
end

% Set auxiliary variables and arrays
t = (0:dt:tmax);   % integration times
nt = numel(t);     % number of times

% Assign a mass to each coordinate
m = repmat(mass,nd,1);

% Integrate trajectory, using velocity-Verlet algorithm
r=zeros(nd,np,nt); v=r; a=r;          % set size of output arrays
r(:,:,1) = r0;                        % positions at t=0
v(:,:,1) = v0;                        % velocities at t=0
a(:,:,1) = force(r0)./m;              % accelerations at t=0
for it = 2:nt
    r(:,:,it) = r(:,:,it-1) + ...
                v(:,:,it-1)*dt + ...
                a(:,:,it-1)*dt^2/2;   % positions at t
    a(:,:,it) = force(r(:,:,it))./m;  % accelerations at t
    v(:,:,it) = (r(:,:,it)-r(:,:,it-1))/dt + ...
                 a(:,:,it)*dt/2;      % velocities at t
end

% Remove singleton dimensions (of size=1)
r = squeeze(r);
v = squeeze(v);
a = squeeze(a);

end % function verlet
