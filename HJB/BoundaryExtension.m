function A = BoundaryExtension(type,B,C,varargin)

    % Linear (constant): ('linear',B,C)
    
    % Linear (variable): ('linear',B,C,[a,b]);
    
    % Quadratic (constant): ('quadratic',B,C,D)
    
    % Quadratic (variable): ('quadratic',B,C,D,[a,b,c])

    switch type
        
       case 'linear'
        
           % Constant spacing
           if nargin == 3
               A  = 2*B - C;
               
           % Variable spacing
           else
               dx = varargin{1};
               A  = (B - C) * (dx(1) / dx(2)) + B;
           end
           
        case 'quadratic'
        
           D = varargin{1};
           % Constant spacing
           if nargin == 4
               A  = 3*B - 3*C + D;
               
           % Variable spacing
           else
               dx = varargin{2};
               a = dx(1);
               b = dx(2);
               c = dx(3);
               A  = ((a + b + c)*((a+b)*c * B - a*(b+c) * C) + a*b*(a+b) * D) / (b*c*(b + c));
           end
                        
    end
        
    
end