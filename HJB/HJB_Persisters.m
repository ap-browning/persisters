function varargout = HJB_Persisters(f,s2,lambda,J,Ufcn,grid,varargin)
%HJB_PERSISTERS Code to solve persister control problem, where the solution
%   V = Vtilde + x.
% This code solves for Vtilde, which is independent of X.
%
%   dX = f{1}(T,U,X,Y,Z) dt + sqrt(s2{1}(T,U,X,Y,Z)) dW
%   dY = f{2}(T,U,X,Y,Z) dt + sqrt(s2{2}(T,U,X,Y,Z)) dW
%   dZ = f{3}(T,U,X,Y,Z) dt + sqrt(s2{3}(T,U,X,Y,Z)) dW
%
% that minimises/maximises the payoff
%
%   J(U) = E [ int_0^T J.C(T,U,X,Y,Z) dt + J.D(XT,YT,ZT) ]
%
% (for fewer than 3 dimensions, Y and/or Z are omitted)
%
% Inputs:
%   f{1}    = @(T,U,X,Y,Z) [used]
%   f{2}    = @(T,U,X,Y,Z)
%   f{3}    = @(T,U,X,Y,Z)
%
%   s2{1}   = @(T,U,X,Y,Z) [unused]
%   s2{2}   = @(T,U,X,Y,Z)
%   s2{3}   = @(T,U,X,Y,Z)
%
%   sxy     = @(T,U,X,Y,Z) [unused]
% 
%   J.C     = @(T,U,Y,Z)
%   J.D     = @(YT,ZT)
%
%   Ufcn    = @(T,Y,Z,V,Vy,Vz,Vyy,Vzz)
%
%
%   grid{1}  = [Tmin, Tmin + dT, ... , Tmax] (equally spaced)
%   grid{2}  = [Ymin, Ymin + dY, ... , Ymax] (variable spaced)
%   grid{3}  = [Zmin, Zmin + dZ, ... , Zmax] (equally spaced)
%
%
% Optional Value/Input pairs:
%
%   Option          (default)   Description
%   --------        --------    --------
%   'SaveEvery'     0           (Must be divisible by size(grid.T))
%   'Bounds'        [0 1]       Bounds
%   'Animate'       false       Show mesh animating solution as solving
%   'Method'        'explicit'  'explicit' or 'heuns'
%
% Example usage:
%   U      = HJB_Persisters(f,s2,lambda,J,Ustar,grid);
%   U      = HJB_Persisters(f,s2,lambda,J,Ustar,grid,'tol',1e-5);
%   [U,V0] = HJB_Persisters(f,s2,lambda,J,Ustar,grid,'tol',1e-5,'maxiters',200);
%   [U,V0] = HJB_Persisters(f,s2,lambda,J,Ustar,grid,'Method','Heuns','maxiters',200);
%   [U,V0] = HJB_Persisters(f,s2,lambda,J,Ustar,grid,'Bounds',[0 0.5]);  
%
% By AP Browning
% Updated 6/11/19

%%%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
addOptional(p,'SaveEvery',  1       );
addOptional(p,'Method',     'explicit');
addOptional(p,'Bounds',     false   );
addOptional(p,'Animate',    false   );
parse(p,varargin{:});

% Load options
SaveEvery   = p.Results.SaveEvery;
Bounds      = p.Results.Bounds;
Animate     = p.Results.Animate;
Method      = p.Results.Method;

if ~Bounds
    Bounded = false;
else
    Bounded = true;
end

tic;

%%%% SETUP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dimension of PDE
dim     = 2;

% Create Ustar
Ustar   = Ufcn;

                           
%%%% GRID / STORAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt          = grid{1}(2) - grid{1}(1);
Nt          = length(grid{1});

dX          = cell([1,dim]);
NX          = zeros(1,dim);
for i = 1:dim
   dX{i}    = diff(grid{i+1}); 
   NX(i)    = length(grid{i+1});
end
NXcell = num2cell(NX);

Xm          = cell([1 dim]);
[Xm{:}]     = ndgrid(grid{2:end});

V           = J.D(Xm{:});

ispoisson   = sum(abs(lambda)) ~= 0;
[~,index_a] = min(abs(grid{3} - -1));   % Index of z = -1
[~,index_b] = min(abs(grid{3} - +1));   % Index of z = +1

% Upwinding scheme (indicator function)
backdiff_y  = false(size(V));

% Used in finite differences for variable spacing
dYprodsum   = (dX{1}(1:end-1) + dX{1}(2:end)) .* dX{1}(1:end-1) .* dX{1}(2:end);
dZprodsum   = (dX{2}(1:end-1) + dX{2}(2:end)) .* dX{2}(1:end-1) .* dX{2}(2:end);

% Create matrix to store (and check size)
if SaveEvery ~= 0
    
    % Continuous problem %%%%%%%%%%
    requestedsizeMB = floor(Nt/SaveEvery)*prod(NX) / (1024^2);
    if requestedsizeMB >= 8000
        error('Requested storage size is greater than 8GB.');
        
    elseif requestedsizeMB > 4000 && requestedsizeMB < 8000
        warning('Requested storage size is between 4GB and 8GB, switching to single precision');
        Ustore  = zeros(floor(Nt/SaveEvery),NXcell{:},'single');
        V       = single(V);
        
    else
        Ustore  = zeros(floor(Nt/SaveEvery),NXcell{:});
        
    end
 
end

%%%% SOLVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = (Nt-1):-1:1
    
    tcur = grid{1}(i+1);
      
    % Solve with specified Method
    switch Method
        case 'explicit'
            V   = V - dt * Fbar(V,tcur);
        case 'heuns'
            dV1  = Fbar(V,tcur);
            V   = V - dt / 2 * (dV1 + Fbar(V - dt * dV1,tcur - dt));
    end
    
    % Save?
    if SaveEvery ~= 0 && mod(i-1,SaveEvery) == 0
        
        Ustore(floor((i-1)/SaveEvery)+1,:,:,:) = GetControl(V,tcur - dt);
        
    end
    
    % Plot, show time
    if Animate
        if mod(i-1,250) == 0
            mesh(grid{2},grid{3},V');
            title(['t = ',num2str(tcur)]);
            drawnow;        
        end
    end


end

if SaveEvery == 0
    Ustore = GetControl(V,tcur - dt);
end

Details.Time = toc;

%%%% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch nargout
    case 1
        varargout{1} = Ustore;
    case 2
        varargout{1} = Ustore;
        varargout{2} = V;
    case 3
        varargout{1} = Ustore;
        varargout{2} = V;
        varargout{3} = Details;
end


%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % PDE (Fbar in Supplementary Material)
    function dV = Fbar(V,t)
        
        % Differentiate
        [Vx,Vxx] = diffV(V);
        
        % Get Control
        U  = Ustar(t,Xm{:},V,Vx{:},Vxx{:});
        
        % Apply bounds
        if Bounded
            U(U < Bounds(1)) = Bounds(1);
            U(U > Bounds(2)) = Bounds(2);
        end
        
        % Payoff
        dV = -J.C(t,U,Xm{:}) - f{1}(t,U,0,Xm{:});
        
        % Update winding indicator
        backdiff_y  = (f{2}(t,U,0,Xm{:})) < 0;
        
        % Infinitesimal generator
        for d = 1:dim
            dV = dV - Vx{d}  .*  f{d+1}(t,U,0,Xm{:}) ...
                    - Vxx{d} .* s2{d+1}(t,U,0,Xm{:}) * 0.5;
        end
        
        % Jump process
        if ispoisson
            % Va(x,y,z) = V(x,y,-1) - V(x,y,z)
            % Vb(x,y,z) = V(x,y,+1) - V(x,y,z)
            Va = -V + V(:,index_a);
            Vb = -V + V(:,index_b);
            dV = dV - lambda(1) * (Va + (Xm{2} + 1) .* Vx{2});
            dV = dV - lambda(2) * (Vb + (Xm{2} - 1) .* Vx{2});
        end
              
    end

    % FINITE DIFFERENCE APPROXIMATION OF DERIVATIVES
    function [Vx,Vxx]  = diffV(V)
        
        Vx      = cell([1,dim]);
        Vxx     = cell([1,dim]);
        
        % Y (use winding for first derivative)
        dVy     = diff(V,1,1);

        dVforward = [dVy ./ dX{1}'; nan(1,NX(2))];
        dVbackward = [nan(1,NX(2)); dVy ./ dX{1}'];
        
        Vx{1}   = (dVbackward .* backdiff_y + dVforward .* (~backdiff_y));

        for j = 1:NX(2)
           
            if isnan(Vx{1}(1,j))
                Vx{1}(1,j) = BoundaryExtension('linear',Vx{1}(2,j),Vx{1}(3,j),dX{1}(1:3));
            end
            
            if isnan(Vx{1}(end,j))
                Vx{1}(end,j) = BoundaryExtension('linear',Vx{1}(end-1,j),Vx{1}(end-2,j),dX{1}(end:-1:end-3));
            end
            
        end
        
        Vxx{1}  = 2 * (-dX{1}(2:end)' .* dVy(1:end-1,:) + dX{1}(1:end-1)' .* dVy(2:end,:)) ./ (dYprodsum');
        Vxx{1}  = [zeros(1,NX(2)); Vxx{1}; zeros(1,NX(2))];
        
        % Z (use central differences + articifial BCs) [not needed for
        % poisson]
        if ~ispoisson
            
            dVz     = diff(V,1,2);

            Vx{2}   = (dX{2}(2:end).^2 .* dVz(:,1:end-1) + dX{2}(1:end-1).^2 .* dVz(:,2:end)) ./ dZprodsum;
            Vx{2}   = [BoundaryExtension('linear',Vx{2}(:,1),Vx{2}(:,2),dX{2}(1:3)), Vx{2}, BoundaryExtension('linear',Vx{2}(:,end),Vx{2}(:,end-1),dX{2}(end:-1:end-3))];
                    
            Vxx{2}  = 2*(-dX{2}(2:end) .* dVz(:,1:end-1) + dX{2}(1:end-1) .* dVz(:,2:end)) ./ dZprodsum;
            Vxx{2}   = [BoundaryExtension('linear',Vxx{2}(:,1),Vxx{2}(:,2),dX{2}(1:3)), Vxx{2}, BoundaryExtension('linear',Vxx{2}(:,end),Vxx{2}(:,end-1),dX{2}(end:-1:end-3))];

        
        else
            
            Vx{2} = zeros(size(V));
            Vxx{2} = zeros(size(V));
        
        end

    end

    % GET CONTROL (FOR STORAGE)
    function U = GetControl(V,t)
       
        [Vx,Vxx] = diffV(V);

        U = Ustar(t,Xm{:},V,Vx{:},Vxx{:});
        if Bounded
            U(U < Bounds(1)) = Bounds(1);
            U(U > Bounds(2)) = Bounds(2);
        end

    end

end