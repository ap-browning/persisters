function varargout = Env_2_Monod(varargin)
% MONOD ENVIRONMENT
%
% USAGE:
%   [f,s2,cov,lambda] = Env_2_Monod()
%       Returns f, s2, scov, lambda for HJB
%
%   [f,S,IC]   = Env_2_Monod()
%       Returns f, S, IC for forward solution where
%           f  = @(t,U,X,Y,Z)
%           S  = @(t,U,X,Y,Z)
%
%   [Z,G]    = Env_2_Monod(T) or Env_1_Monod(T,seed) or Env_1_Monod(T,true)
%       Returns generated environment, Z, and growth rate, G
%
%   [G]      = Env_2_Monod(Z,'mu')
%       Maps input environment Z to equivalent growth rate, G
%
%   [Z]      = Env_2_Monod(G,'inv')
%       Maps input G to equivalent monod environment, Z
%

% General Parameters
addpath('../');
Parameters;

% Specific Parameters
mu_max      = 8;
delta       = 2;

Ks          = 1;
theta       = 0.1;
kappa       = 0.3;

% Cutoff value
Zmin    = -0.5;

% Noise to Growth function (& Inverse)
mu      = @(Z) mu_max .* max(Zmin,Z) ./ (Ks + max(Zmin,Z)) - delta;
muinv   = @(G) (-delta * Ks - G * Ks) ./ (delta + G - mu_max);

% Environment functions
fEnv    = @(T,Z) theta * (Ks - Z);
vEnv    = @(T,Z) kappa^2;

% Z IC
z0      = Ks;

% Z limits
EnvZLims    = [-0.5,2.5];


if nargin == 0
%% Build functions for HJB or HJB_Forward

    % Load state equations
    StateEquations;

    f{3}    = @(T,U,X,Y,Z) fEnv(T,Z);

    s2{3}   = @(T,U,X,Y,Z) vEnv(T,Z); 
      
    lambda  = [0,0];
                        
    switch nargout
        case 3
            % For HJB_Forward
            varargout{1} = F;
            varargout{2} = S;
            varargout{3} = IC;
    
        case 4
        % For HJB
            varargout{1} = f;
            varargout{2} = s2;
            varargout{3} = sxy;
            varargout{4} = lambda;
            
        case 5
        % For HJB (suggested z limits)
            varargout{1} = f;
            varargout{2} = s2;
            varargout{3} = sxy;
            varargout{4} = lambda;
            varargout{5} = EnvZLims;
            
    end
    
else
    if nargin == 1 || ~(strcmp(varargin{2},'inv') || strcmp(varargin{2},'mu'))
        
        % Set seed
        if nargin == 2
            if varargin{2} == true
                rng(EnvironmentSeed('monod'));
            else
                rng(varargin{2});
            end
        else
            rng(EnvironmentSeed);
        end
        
        T = varargin{1};
        dt = T(2) - T(1);
        Z = zeros(size(T));
        Z(1) = z0;
        
        for i = 2:length(T)
            
            Z(i) = Z(i-1) + dt * fEnv(T(i-1),Z(i-1)) + sqrt(dt) * sqrt(vEnv(T(i-1),Z(i-1))) * randn(1,1);
            
        end
        
        G = mu(Z);
        
        varargout{1} = Z;
        varargout{2} = G;
        
    else
%% Mu and InverseMu function
    
        if strcmp(varargin{2},'inv')
            varargout{1} = muinv(varargin{1});
        elseif strcmp(varargin{2},'mu')
            varargout{1} = mu(varargin{1});
        end
        
    end
    
end

% Shuffle generator
rng('shuffle');