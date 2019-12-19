%% STATE EQUATIONS

% dX = F(T,U,X) dt + S(T,U,X) dW

% dn  = f_{1} dt + S_{11} dW1 + S_{12} dW2
% dth = f_{2} dt + S_{21} dW1 + S_{22} dW2

%%%% FOR HJB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% f{i}  = F_{i}
f{1}    = @(T,U,X,Y,Z) (1 - (1 - eps) * Y) .* mu(Z) - 0.5 * sigma^2 * (1 - Y).^2 - 0.5  * eta^2 * Y.^2;
f{2}    = @(T,U,X,Y,Z) (umin + U) .* (1 - Y) - v * Y + (-1 + Y) .* Y .* (eta^2 * Y + mu(Z) * (1 - eps) + (-1 + Y) * sigma^2);

% s2{i} = (SS')_{ii}
s2{1}    = @(T,U,X,Y,Z) sigma^2 * (1 - Y).^2 + eta^2 * Y.^2;
s2{2}    = @(T,U,X,Y,Z) sigma^2 * (1 - Y).^2 .* Y.^2 + eta^2 * (1 - Y).^2 .* Y.^2;

% sxy   = (SS')_{12} = (SS')_{21}
sxy     = @(T,U,X,Y,Z) -sigma^2 * (1 - Y).^2 .* Y + eta^2 * (1 - Y) .* Y.^2;

%%%% FOR SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F       = @(T,U,X,Y,Z) [f{1}(T,U,X,Y,Z); f{2}(T,U,X,Y,Z)];
    
S       = @(T,U,X,Y,Z) [ (1 - Y) * sigma,      eta * Y;
                        -(1 - Y) .* Y * sigma, eta * (1 - Y) .* Y];