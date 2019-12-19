function varargout = Env_1_Constant(varargin)
% MONOD ENVIRONMENT
%
% USAGE:
%   [f,s2,cov,lambda] = Env_1_Constant()
%       Returns f, s2, scov, lambda for HJB
%
%   [f,S,IC]   = Env_1_Constant()
%       Returns f, S, IC for forward solution where
%           f  = @(t,U,X,Y,Z)
%           S  = @(t,U,X,Y,Z)
%
%   [Z,G]    = Env_1_Constant(T) or Env_1_Constant(T,seed) or Env_1_Constant(T,true)
%       Returns generated environment, Z, and growth rate, G
%
%   [G]      = Env_1_Constant(Z,'mu')
%       Maps input environment Z to equivalent growth rate, G
%
%   [Z]      = Env_1_Constant(G,'inv')
%       Maps input G to equivalent monod environment, Z
%

% General Parameters
addpath('../');
Parameters;

% Specific Parameters
muG         = 2.0;

% Noise to Growth function (& Inverse)
mu          = @(Z) Z;
muinv       = @(G) G;

% Environment functions
fEnv        = @(T,Z) 0;
vEnv        = @(T,Z) 0;

% Z IC
z0          = muG;

% Z limits
EnvZLims    = [-3,3];

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
            if varargin{2}
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