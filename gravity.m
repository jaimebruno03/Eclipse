function f = gravity( m, r )

% Finds the gravitational force between several masses
% J.M.Soler, Jan.2011
%
% Input:
%   m(np)    : mass of the np particles, in kg
%   r(nd,np) : position of np particles in nd dimensions, in m
% Output:
%   f(nd,np) : force on the np particles in each of the nd dimensions, in N

% Set gravitational constant
G = 6.67428e-11;  % in N*m^2/kg^2

% Find number of particles and space dimensions
[nd,np] = size(r);

% Find forces
f = zeros(nd,np);
for ip = 1:np-1
    for jp = ip+1:np
        rij = r(:,jp) - r(:,ip);           % vetor from ip to jp
        dij = norm(rij);                   % distance between ip and jp
        uij = rij/dij;                     % unit vector from ip to jp
        fij = -G*m(ip)*m(jp)/dij^2 * uij;  % force of ip on jp
        f(:,jp) = f(:,jp) + fij;
        f(:,ip) = f(:,ip) - fij;
    end
end

end % function gravity
