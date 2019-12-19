function grid = HJB_CreateGrid(variable,varargin)
%HJB_CreateGrid Creates grid for HJB equations
%
% INPUTS:
%   variable    'y' or 'z'
%   
% ADDITIONAL INPUTS (variable = 'y')
%   Nlog        number in log space component
%   Nlin        number in linear space component
%   lmin        
%   lmax
%
% grid is 
%        [0  lmin    lmax   1-lmax   1-lmin  1];
%    Ex: [0  1e-7     0.1     0.9    1-1e-7  1]
%              |  log  |  lin  |  log   |
%
% ADDITIONAL INPUTS (variable = 'z')
%   N           number total (linear spacing)
%   lims        limits

switch variable
    case 'y'
        
        Nlog        = varargin{1};
        Nlin        = varargin{2};
        lmin        = varargin{3};
        lmax        = varargin{4};
        
        Ylogvec     = logspace(log10(lmin),log10(lmax),Nlog);
        dYlogvecend = Ylogvec(end) - Ylogvec(end-1);
        Ylinvec     = linspace(lmax + dYlogvecend,1 - lmax - dYlogvecend,Nlin);

        grid        = [0, Ylogvec, Ylinvec, 1 - Ylogvec(end:-1:1), 1];         % Y
        
    case 'z'

        N           = varargin{1};
        lims        = varargin{2};
        
        grid        = linspace(lims(1),lims(2),N);

end