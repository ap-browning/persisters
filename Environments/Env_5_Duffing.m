function varargout = Env_5_Duffing(varargin)
% MONOD ENVIRONMENT
%
% USAGE:
%   [f,s2,cov,lambda] = Env_5_Duffing()
%       Returns f, s2, scov, lambda for HJB
%
%   [f,S,IC]   = Env_5_Duffing()
%       Returns f, S, IC for forward solution where
%           f  = @(t,U,X,Y,Z)
%           S  = @(t,U,X,Y,Z)
%
%   [Z,G]    = Env_5_Duffing(T) or Env_5_Duffing(T,seed) or Env_5_Duffing(T,true)
%       Returns generated environment, Z, and growth rate, G
%
%   [G]      = Env_5_Duffing(Z,'mu')
%       Maps input environment Z to equivalent growth rate, G
%
%   [Z]      = Env_5_Duffing(G,'inv')
%       Maps input G to equivalent monod environment, Z
%

% General Parameters
addpath('../');
Parameters;

% Specific Parameters
muG         = 2;
muS         = -2;

a           = 20;
b           = 0.4;
nu          = 0.5;

% Noise to Growth function (& Inverse)
mu          = @(Z) muS + (muG - muS) * Z;
muinv       = @(G) (muS - G) / (muS - muG);

% Environment functions
fEnv    = @(T,Z) a * Z .* (1 - Z) .* (Z - b);
vEnv    = @(T,Z) nu^2;

% Z IC
z0      = 1;

% Z limits
EnvZLims    = [-0.5,1.5];

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
                rng(EnvironmentSeed('duffing'));
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