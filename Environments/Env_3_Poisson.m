function varargout = Env_3_Poisson(varargin)
% MONOD ENVIRONMENT
%
% USAGE:
%   [f,s2,cov,lambda] = Env_3_Poisson()
%       Returns f, s2, scov, lambda for HJB
%
%   [f,S,IC]   = Env_3_Poisson()
%       Returns f, S, IC for forward solution where
%           f  = @(t,U,X,Y,Z)
%           S  = @(t,U,X,Y,Z)
%
%   [Z,G]    = Env_3_Poisson(T) or Env_1_Monod(T,seed) or Env_1_Monod(T,true)
%       Returns generated environment, Z, and growth rate, G
%
%   [G]      = Env_3_Poisson(Z,'mu')
%       Maps input environment Z to equivalent growth rate, G
%
%   [Z]      = Env_3_Poisson(G,'inv')
%       Maps input G to equivalent monod environment, Z
%

% General Parameters
addpath('../');
Parameters;

% Specific Parameters
tau         = [9.5,0.5];
lambda      = 1 ./ tau;

% Noise to Growth function (& Inverse)
mu          = @(Z) 2   * Z;
muinv       = @(G) 0.5 * G;

% Environment functions
fEnv    = @(T,Z) 0 * Z;
vEnv    = @(T,Z) 0 * Z;

% Z IC
z0      = 1;

% Z limits
EnvZLims = [-1,1];

if nargin == 0
%% Build functions for HJB or HJB_Forward

    % Load state equations
    StateEquations;

    f{3}    = @(T,U,X,Y,Z) fEnv(T,Z);

    s2{3}   = @(T,U,X,Y,Z) vEnv(T,Z); 
    
                        
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
                rng(EnvironmentSeed('poisson'));
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
            
            if Z(i-1) > 0
                dP  = (poissrnd(lambda(1) * dt) > 0) * -2;
            else
                dP  = (poissrnd(lambda(2) * dt) > 0) * +2;
            end    
            
            Z(i) = Z(i-1) + dt * fEnv(T(i-1),Z(i-1)) + sqrt(dt) * sqrt(vEnv(T(i-1),Z(i-1))) * randn(1,1) + dP;

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