function [sys,x0,str,ts,simStateCompliance] = RP_chrom_try_2(t,x,u,flag,gridsize)

switch flag,

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0,
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(gridsize);
  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1,
    sys=mdlDerivatives(t,x,u,gridsize);
  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2,
    sys=mdlUpdate(t,x,u);
  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3,
    sys=mdlOutputs(t,x,u,gridsize);
  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4,
    sys=mdlGetTimeOfNextVarHit(t,x,u);
  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9,
    sys=mdlTerminate(t,x,u);
  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(gridsize)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%

par = struct(...
    'components',3,...          % Number of components
    'nconc',3);                 % Number of concentrations (mobile phase concentration, intrapore concentration, concentration of adsorbed specie)


sizes = simsizes;

sizes.NumContStates  = par.nconc*par.components*gridsize;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = par.nconc*par.components*gridsize;
sizes.NumInputs      = par.components;
sizes.DirFeedthrough = 1;   % to understand whether inputs affect outputs (1=yes, 0=no)
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%

x0  = zeros(par.nconc*par.components*gridsize ,1); 
% x0 = [] 
%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,gridsize)

%% parameters
par = struct(...  
    'components',3,...          % Number of components
    'length', 500E-2,...            % Length of column [dm]
    'velocity', 20/1000/0.03*10,...        % Linear velocity [dm/h] Q=20 L/h A=0.03 m^2 velocity=Q/A
    'dax', 0.0113,...             % Axial dispersion [dm^2/h] 
    'epsbed',  0.35,...         % Bed porosity
    'epspar',  0.8,...          % Particle porosity
    'keff',  [0.4172 0.5606 0.3443],...  % Effective film transfer [dm/h] 
    'rad', 10E-6,...            % Particle radius [dm]
    'ads', [1E-07*0.9,1E-07*0.9,6*1.1E2]*3600E-6,...    % Langmuir adsorption coefficient [dm^3/g/h]
    'des', [4E09*1.1,4E09*1.1,1*0.9E-3]*3600,...        % Langmuir desorption coefficient [1/h]
    'qmax', [1E-6*0.9,1E-6*0.9,70*1.1]);      % Langmuir maximum concentration %¤[g/dm^3]
            
%% parameters ok
% par = struct(...  
%    'components',3,...          % Number of components
%    'length', 500E-2,...            % Length of column [dm]
%    'velocity', 0.06*3600E-2,...        % Pump velocity [dm/h]
%    'dax', 0.0113,...             % Axial dispersion [dm^2/h] 
%    'epsbed',  0.35,...         % Column porosity
%    'epspar',  0.8,...          % Particle porosity
%    'keff',  [0.4172 0.5606 0.3443],...  % Effective film transfer [dm/h] 
%    'rad', 0.045E-2,...            % Particle radius [dm]
%    'ads', [1E-07,1E-07,5E2]*3600E-6,...           % Langmuir adsorption coefficient [dm^3/g/h]
%    'des', [4E09,4E09,2E-3]*3600,...            % Langmuir desorption coefficient [1/h]
%    'qmax', [1E-6,1E-6,60]);      % Langmuir maximum concentration %¤[g/dm^3]
%% implementation finite elements
nconc = 3; %number of concentrations
h = par.length/(gridsize-1); % Length of each cell
%mass matrix
M =   h/6*diag(ones(gridsize-1,1),1) +... % Upper diagonal
    2*h/3*diag(ones(gridsize,1)) +...     % Main diagonal
      h/6*diag(ones(gridsize-1,1),-1);    % Lower diagnal
M(gridsize,gridsize)=h/3;                      % Correct last entry
M(1,1)=1;                            % Prepare Dirichlet BC
M(1,2)=0;
M=sparse(M); %[dm]
%convection matrix
C = 1/2*diag(ones(gridsize-1,1), 1) + ... % Upper diagonal
      0*diag(ones(gridsize  ,1))    + ... % Main diagonal
   -1/2*diag(ones(gridsize-1,1),-1);      % Lower diagnal
C(1,1)     = -1/2;                   % Correct first entry
C(gridsize,gridsize) =  1/2;                   % Correct last entry
C=sparse(C);
%
A = -1/h*diag(ones(gridsize-1,1), 1) + ...% Upper diagonal
     2/h*diag(ones(gridsize  ,1))    + ...% Main diagonal
    -1/h*diag(ones(gridsize-1,1),-1);     % Lower diagnal
A(1,1)    =1/h;                      % Correct first entry
A(gridsize,gridsize)=1/h;                      % Correct last entry
A=sparse(A); %[1/dm]
%% states

for k=1:nconc
    for  i=1:par.components
        for j=1:gridsize
            c(j,i,k) = x(gridsize*par.components*(k-1)+gridsize*(i-1)+j);
        end
    end
end

%% inputs

for i=1:par.components
    conc(i) = u(i);
end

%% equations

    for j=1:par.components
    
        % M*dc/dt = -(vC+dA)c -(1-epsbed)/epsbed*keff*3/rad*M*(c-c_p) +BC
        dcdt(:,j,1) = -(par.velocity/par.epsbed*C + par.dax*A)*...
            c(:,j,1) - (1-par.epsbed)/par.epsbed*par.keff(j)*3.0/...
            par.rad*M*(c(:,j,1)-c(:,j,2));

        % Boundary condition
        dcdt(1,j,1) = dcdt(1,j,1) - ...
            par.velocity/par.epsbed*(c(1,j,1)-conc(j));

    	dcdt(:,j,1) = M\dcdt(:,j,1);         % Incorporate M
        
        % Langmuir model
        sum=ones(gridsize,1);
        for k=1:par.components
            sum=sum-c(:,k,3)/par.qmax(k);
        end
        dcdt(:,j,3) = par.ads(j)*par.qmax(j)*sum.*c(:,j,2) - ...
            par.des(j)*c(:,j,3);
        
        % Pore concentration
        dcdt(:,j,2) = -(1-par.epspar)/par.epspar*dcdt(:,j,3) + ...
            par.keff(j)*3/par.rad/par.epspar*(c(:,j,1)-c(:,j,2));     

    end

dcdt = reshape(dcdt,[],1,1); %dcdt as a column vector
sys = [dcdt(:)];

% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];


% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%

function sys=mdlOutputs(t,x,u,gridsize)
ncomp = 3;
nconc = 3;
% x=reshape(x,[],1);
sys = eye(ncomp*gridsize*nconc)*x;


% function sys=mdlOutputs(t,x,u)
% % dim = 101; 
% % x=reshape(x,[],1);
% sys = [x(:,1); x(:,2)];

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate