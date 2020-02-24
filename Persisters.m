%% VARIABLE PERSISTER PRODUCTION
%  Used to generate Figs 2, 3, S1, S2, S4 and S5

%% LOAD CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ENVIRONMENTS
Envs = {'Env_1_Constant'            , ...
        'Env_2_Monod'               , ...
        'Env_3_Poisson'             , ...
        'Env_4_OrnsteinUhlenbeck'   , ...
        'Env_5_Duffing'            };

% ADD REQUIRED FILES TO PATH
  addpath('Environments');
  addpath('HJB');
    
  
%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ? ENVIRONMENT ?
% Environment to solve (Poisson is much faster)
  Env     = Envs{5};

% ? LOAD RESULTS ? 
% Option to save/load results stored in "Environments/%Env%.mat"
%   If 'Load' is true, the HJB equation will not solve, but rather load
%   solution from the last time 'Save' was turned on for the respective
%   environment.
  Save    = false;
  Load    = false;

% ? PROBLEM PARAMETERS ?
% Final time
  tf      = 10;
% Control weighting
  alpha   = -100;
% Control bounds
  bounds  = [0,0.1];

  
%% SOLVE HJB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LOAD ENVIRONMENT EQUATIONS
  [f,s2,cov,lambda,EnvZLims] = feval(Env);

% PAYOFF FUNCTION
  J.C       = @(T,U,Y,Z) alpha * U.^2;
  J.D       = @(YT,ZT)   0 * YT;
  Ufcn      = @(T,Y,Z,V,Vy,Vz,Vyy,Vzz) -Vy .* (1 - Y) / (2 * alpha);

% SPATIAL DISCRETIZATION
  Nlog      = 49;
  Nlin      = 30;
  lmin      = 1e-4;
  lmax      = 0.1;
  Nz        = 100;
  if strcmp(Env,'Env_3_Poisson')
      Nz  = 2;
  elseif strcmp(Env,'Env_4_OrnsteinUlenbeck')
      Nz  = 200;
  end

  % Time (s)
  grid{1}     = linspace(0,tf,300000);
  
  % y
  grid{2}     = HJB_CreateGrid('y',Nlog,Nlin,lmin,lmax);
  
  % z
  grid{3}     = HJB_CreateGrid('z',Nz,EnvZLims);
  
% SAVE EVERY
%  Generally, the timestep for the SDE forward problem will be much larger
%  than that for solving the PDE. Only save every 'SaveEvery' time points
%  in the PDE, starting at time s = 0 in the PDE solution.
  SaveEvery = 60;

% LOAD ?
if Load
    load(['Environments/',Env,'.mat']);
    
% SOLVE HJB
else
    [U,V0,HJBspec] = HJB_Persisters(f,s2,lambda,J,Ufcn,grid, ...
            'SaveEvery',SaveEvery   , ...
            'Bounds',bounds        );
end

% SAVE ?
if Save
    save(['Environments/',Env,'.mat'],'Env','U','V0','SaveEvery','grid','-v7.3');
end

% VISUALIZE SOLUTION
    figure(1); clf; set(0,'defaultAxesFontSize',12)

    % V0
    subplot(1,2,1);
        surf(grid{2},grid{3},V0');
        xlabel('y'); ylabel('z'); zlabel('V');
        title('V');

    % Phi_0
    subplot(1,2,2);
        surf(grid{2},feval(Env,grid{3},'mu'),reshape(U(1,:,:),length(grid{2}),length(grid{3}))');
        xlabel('\theta_0'); ylabel('\mu_0'); zlabel('\phi_0');
        title('Variable persister rate');

    
%% SOLVE SDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TIME GRID FOR SDE
  T           = [grid{1}(1:SaveEvery:end),tf];

% LOAD ENVIRONMENT EQUATIONS
  [F,S,IC]    = feval(Env);

% SIMULATE INDEPENDENT ENVIRONMENT SDE (GENERATE ENVIRONMENT)
  [Z,G]       = feval(Env,T,true);
%                            ^^ 'true' to use seed from main document.
%                               Remove for random, or set to another seed.

% SIMULATE UNDER ANTIBIOTICS
%  G           = -2 * ones(size(T));
%  Z           = feval(Env,G,'inv');

% SIMULATE STATE EQUATIONS, USING ENVIRONMENT AND CONTROL
  [Y,Uy]      = HJB_Forward_Persisters(T,Z,F,S,IC,U,grid,true);
%                'true' to use seed from main document ---^^
%                Remove for random, or set to another seed.

% VISUALIZE SOLUTION
    figure(2); clf; set(0,'defaultAxesFontSize',12);
                    set(0,'DefaultLineLineWidth', 2);
        
    % Environment
    subplot(1,3,1);
        plot(T,G);
        axis([0,10,-3,3]);
        xlabel('t'); ylabel('\mu_t');
        title('Environment');
        
    % Control
    subplot(1,3,2);
        plot(T,Uy);
        axis([0,10,-0.01,0.11]);
        xlabel('t'); ylabel('\phi_t');
        title('Variable persister production');
        
    % Populations
    subplot(1,3,3);
        yyaxis left;
        semilogy(T,exp(Y(1,:)));
        axis([0,10,10.^[-1,8]]);
        ylabel('n_t');
        
        yyaxis right;
        semilogy(T,Y(2,:),'LineWidth',2);
        axis([0,10,10.^[-7,0]]);
        ylabel('\theta_t');
        
        xlabel('t');
        title('Variable persister production');
  
