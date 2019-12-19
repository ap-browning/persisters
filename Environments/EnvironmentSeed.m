function seed = EnvironmentSeed(varargin)
%ENVIRONMENTSEED generate seed for environment generation
%
%Usage:
%   seed = EnvironmentSeed()  generates random seed
%   seed = EnvironmentSeed(1) returns seed for standard environment 1
%   seed = EnvironmentSeed('gaussian') is equivalent to EnvironmentSeed(1)

rng('shuffle');

% Stored seeds
EnvironmentStored = [ 863914;       % Ornstein
                      179054;       % Monod
                      593766;       % Poisson
                      599356  ];    % Duffing

switch nargin

    case 0
        seed = randi(1e6);
        
    case 1
        env  = varargin{1};
        switch env
            case 'ornstein'
                env = 1;
            case 'monod'
                env = 2;
            case 'poisson'
                env = 3;
            case 'duffing'
                env = 4;
        end
        seed = EnvironmentStored(env);
    
end