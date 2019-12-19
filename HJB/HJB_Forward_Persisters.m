function varargout = HJB_Forward_Persisters(T,Z,F,S,IC,Ustore,grid,varargin)

% INPUTS:
%   T       = [Tmin, Tmin + dT, ... , Tmax]
%   Z       = [z_tmin, z_(tmin + dT), ....]
%   F       = @(T,U,X,Y,Z) [f^(ntilde)(T,U,X,Y,Z); f^(theta)(T,U,X,Y,Z)];
%   S       = @(T,U,X,Y,Z) (2 * 2 volatility scaling matrix)
%   IC      = [ntilde0; theta0];
%   Ustore  = U output from HJB code
%   grid    = grid used in HJB code
%   seed (option) (2 * 1 vec) set to 'true' to use seeds in main paper, else specify or
%       leave out.
%
% OUTPUTS:
%   Y               = [ntilde_Tmin, ntilde_(Tmin + dT),...; 
%                      theta_Tmin,  theta_(Tmin + dT), ...]; 
%   U (optional)    = [phi_Tmin, phi_(Tmin + dT), ...];

%%%% SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y       = zeros(2,length(T));
y       = IC(:);
Y(:,1)  = y;

dt      = T(2) - T(1);

% Generate brownian motion - todo : take seed as input
switch nargin
    
    % No seed
    case 7
        B(1,:)  = randn(1,length(T));
        B(2,:)  = randn(1,length(T));
        
    case 8
        if varargin{1}
            rng(365867);
            B(1,:)  = randn(1,length(T));
            rng(741941);
            B(2,:)  = randn(1,length(T));
        else
            rng(varargin{1}(1));
            B(1,:)  = randn(1,length(T));
            rng(varargin{1}(2));
            B(2,:)  = randn(1,length(T));
        end
        
end

rng('shuffle');

gridY   = grid{2};
gridZ   = grid{3};
[~,NY,NZ] = size(Ustore);

U       = zeros(size(T));

%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 2:length(T)

    % Get current environment
    z           = Z(i-1);
    
    % Get current control
    u           = interpn(gridY,gridZ,reshape(Ustore(i-1,:,:),NY,NZ),y(2),z);

    Fy          = F(T(i-1),u,y(1),y(2),z);
    Sy          = S(T(i-1),u,y(1),y(2),z);

    y           = y +  dt * Fy + sqrt(dt) * Sy * B(:,i-1);

    % Reflecting boundary for theta
    y(2)        = abs(y(2));

    Y(:,i)    = y;
    U(i-1)    = u;
    
end

switch nargout
    
    case 1
        varargout{1} = Y;
        
    case 2
        varargout{1} = Y;
        varargout{2} = U;
        
end
    

end