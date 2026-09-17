function main_directCollocation_Absolute_DoublePendulum_OpenSim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Import CasADi and definition of Implicit or Explicit formulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path to where CasADi folder is located
%path('root_to_folder');
import casadi.*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Collocation points and definition of constant matrices %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Degree of interpolating polynomial
d = 3;

% Get collocation points
tau_root = collocation_points(d, 'legendre');
tau_root=[0 tau_root];

% Coefficients of the collocation equation
C = zeros(d+1,d+1);

% Coefficients of the continuity equation
D = zeros(d+1, 1);

% Coefficients of the quadrature function
B = zeros(d+1, 1);

% Construct polynomial basis
for j=1:d+1
  % Construct Lagrange polynomials to get the polynomial basis at the collocation point
  coeff = 1;
  for r=1:d+1
    if r ~= j
      coeff = conv(coeff, [1, -tau_root(r)]);
      coeff = coeff / (tau_root(j)-tau_root(r));
    end
  end
  % Evaluate the polynomial at the final time to get the coefficients of the continuity equation
  D(j) = polyval(coeff, 1.0);

  % Evaluate the time derivative of the polynomial at all collocation points to get the coefficients of the continuity equation
  pder = polyder(coeff);
  for r=1:d+1
    C(j,r) = polyval(pder, tau_root(r));
  end

  % Evaluate the integral of the polynomial to get the coefficients of the quadrature function
  pint = polyint(coeff);
  B(j) = polyval(pint, 1.0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Set constant parameter values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time horizon
T = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of coordinates and torque controls
nq=6;
ncT=6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Geometric and dynamics parameters
fixedPoints = [
    360, 0, 0;   % O1
    180, 0, 311.77;   % O2
    -180, 0, 311.77;   % O3
    -360, 0, 0;   % O4
    -180, 0, -311.77;   % O5
    180, 0, -311.77;    % O6
];


 auxdata.g=9.81; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Set model variables and model equations  %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Declare model variables
q1 = MX.sym('q1');
q2 = MX.sym('q2');
q3 = MX.sym('q3');
q4 = MX.sym('q4');
q5 = MX.sym('q5');
q6 = MX.sym('q6');
q1dot = MX.sym('q1dot');
q2dot = MX.sym('q2dot');
q3dot = MX.sym('q3dot');
q4dot = MX.sym('q4dot');
q5dot = MX.sym('q5dot');
q6dot = MX.sym('q6dot');
x=[q1 q1dot q2 q2dot q3 q3dot q4 q4dot q5 q5dot q6 q6dot];
uT1 = MX.sym('uT1');
uT2 = MX.sym('uT2');
uT3 = MX.sym('uT3');
uT4 = MX.sym('uT4');
uT5 = MX.sym('uT5');
uT6 = MX.sym('uT6');
uT=[uT1 uT2 uT3 uT4 uT5 uT6];
ua_q1=MX.sym('ua_q1');
ua_q2=MX.sym('ua_q2');
ua_q3=MX.sym('ua_q3');
ua_q4=MX.sym('ua_q4');
ua_q5=MX.sym('ua_q5');
ua_q6=MX.sym('ua_q6');
ua=[ua_q1 ua_q2 ua_q3 ua_q4 ua_q5 ua_q6];

% Model equations
F=external('F','stewardplatform.dll');
xdot = [q1dot;ua_q1; q2dot; ua_q2; q3dot; ua_q3;q4dot;ua_q4; q5dot; ua_q5; q6dot; ua_q6];

s.uT=1e6;
s.ua=100;

% Cost function
%L = sum((uT).^2)+0.01*sum((ua).^2);
% L = sum(ua.^2)+0.01*sum((uT).^2);
L = sum((uT).^2)+sum(ua.^2);


% L = sum(uT.*(x(8:2:nq*2)))+0.01*sum((ua).^2); 


% Boundary conditions
%fixed initial states
initial_states_lb=[        0;   0;           0;   0;     0;   0; 0;    0;   1;   0; 0; 0]; %q1 q1dot q2 q2dot q3 q3dot t4 t4dot t5 t5dot t6 t6dot
initial_states_ub=[        0;   0;           0;   0;     0;   0; 0;    0;   1;   0; 0; 0];
initial_states_ig=[        0;   0;           0;   0;     0;   0; 0;    0;   1;   0; 0; 0];
%fixed final states
final_states_lb=[-14*2*pi/360;   0; 10*2*pi/360;   0; -10*pi/180;   0; 0;    0;   1;   0; 0; 0];
final_states_ub=[-14*2*pi/360;   0; 10*2*pi/360;   0; -10*pi/180;   0; 0;    0;   1;   0; 0; 0];
final_states_ig=[-14*2*pi/360;   0; 10*2*pi/360;   0; -10*pi/180;   0; 0;    0;   1;   0; 0; 0];
%state bounds during the movement
states_lb=[              -2*pi; -30;         -2*pi; -30;   -2*pi; -30; -5; -30; 0; -30; -5; -30]; 
states_ub=[               2*pi;  30;          2*pi;  30;    2*pi;  30;  5;  30; 5;  30;  5;  30];
states_ig=[                0;   0;           0;   0;     0;   0;  0;   0;   1;   0;  0;  0];
%control bounds
controls_lb=-1e8*ones(nq,1)/s.uT;
controls_ub= 1e8*ones(nq,1)/s.uT;
controls_ig=zeros(nq,1);
%acceleration bounds (implicit form)
accelerations_lb=-100*ones(nq,1)/s.ua;
accelerations_ub= 100*ones(nq,1)/s.ua;
accelerations_ig= zeros(nq,1);

% Control discretization
N = 20; % number of control intervals (20 intervals/second reasonable)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Continuous time dynamics
f = Function('f', {x, uT, ua}, {xdot, L});

h = T/N;

%% time grid
tgrid = linspace(0, T, N+1); % tgrid only contains mesh points
for i=1:4
    dtime(i)=tau_root(i)*T/N;
end
for i=1:N
    tgrid_ext([((i-1)*4+1):1:i*4])=[tgrid(i)+dtime];
end
tgrid_ext(end+1)=T; %tgrid_ext contains mesh points and collocation points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Definition and solve of NLP %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% "Lift" initial conditions
X0 = MX.sym('X0',nq*2);
w = {w{:}, X0};

lbw = [lbw; initial_states_lb];      %lower bounds of the initial state
ubw = [ubw; initial_states_ub];     %upper bounds of the initial state
w0 =  [ w0; initial_states_ig];      %initial guess of initial state


%% Define collocation points for the states at the first interval, their bounds and initial conditions
Xkm1j={};
for j=1:d
    Xkm1j{j} = MX.sym(['X_0_' num2str(j)], nq*2);
    w = {w{:}, Xkm1j{j}};
    lbw = [lbw; states_lb];         %lower bound for x_kj (j collocation point of x_k-1)    %%%%%%%%%%%%%%%%%%%%%%%%%
    ubw = [ubw; states_ub];         %upper bound for x_kj  (j collocation point of x_k-1)   %%%%%%%%%%%%%%%%%%%%%%%%%
    w0 = [w0; states_ig];   
end

%% Define controls at the first point
Ukm1= MX.sym(['U_0'],ncT);   %torque controls
w = {w{:}, Ukm1};   
lbw = [lbw; controls_lb];       %lower bound for u_k
ubw = [ubw; controls_ub];       %upper bound for u_k
w0 = [w0; controls_ig];

Uakm1 = MX.sym(['Ua_0'],nq);    %acceleration controls
w = {w{:}, Uakm1};    
lbw = [lbw; accelerations_lb];       %lower bound for u_k
ubw = [ubw; accelerations_ub];       %upper bound for u_k
w0 =  [w0;  accelerations_ig];

% Formulate the NLP
Xkm1 = X0;
for k=1:N-1
    % New NLP variable for states
    Xk = MX.sym(['X_' num2str(k)], nq*2);
    w = {w{:}, Xk};
    lbw = [lbw;  states_lb];             
    ubw = [ubw;  states_ub];              
    w0 = [w0; states_ig];
    
    % New NLP variable for the torque control
    Uk = MX.sym(['U_' num2str(k)],ncT);
    w = {w{:}, Uk};   
    lbw = [lbw; controls_lb];        %lower bound for u_k
    ubw = [ubw; controls_ub];        %upper bound for u_k
    w0 = [w0; controls_ig];          %initial guess for u_k

% New NLP variable for the acceleration control
Uak = MX.sym(['U_a' num2str(k)],nq);
w = {w{:},Uak};
lbw = [lbw; accelerations_lb];     %lower bound for ua_k
ubw = [ubw; accelerations_ub];      %upper bound for ua_k
w0 =  [w0;  accelerations_ig];        %initial guess for ua_k

% State and control at collocation points (control is not a design variable at the collocation points)
    Xkj = {};
    Ukj={};
    Uakj={};
    for j=1:d
        Xkj{j} = MX.sym(['X_' num2str(k) '_' num2str(j)], nq*2);
        w = {w{:}, Xkj{j}};
        lbw = [lbw; states_lb];         %lower bound for x_kj (j collocation point of x_k)   %%%%%%%%%%%%%%%%%%%%%%%%%
        ubw = [ubw; states_ub];         %upper bound for x_kj  (j collocation point of x_k)   %%%%%%%%%%%%%%%%%%%%%%%%%
        w0 = [w0; states_ig];             %initial guess for x_kj (j collocation point of x_k)  %%%%%%%%%%%%%%%%%%%%%%%%%
       
        %Interpolate controls with lagrange polynomials at the collocation points
        Ukj{j}=MX(ncT,1);
        Uakj{j}=MX(nq,1);
      
        Ukj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Ukm1 Uk]);
        Uakj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Uakm1 Uak]);
    end
        
    % Loop over collocation points
    Xk_end = D(1)*Xkm1;
    for j=1:d
       % Expression for the state derivative at the collocation point
       xp = C(1,j+1)*Xkm1;
       for r=1:d
           xp = xp + C(r+1,j+1)*Xkm1j{r};
       end
      
       % Append collocation equations

        [fj qj]=f(Xkm1j{j},Ukj{j},Uakj{j}*s.ua);
        g = {g{:}, h*fj - xp};               %build defect constraints
        lbg = [lbg; zeros(2*nq,1)];             %lower bounds for defect constraints
        ubg = [ubg; zeros(2*nq,1)];             %upper bounds for defect constraints
      
       % Add contribution to the end state
       Xk_end = Xk_end + D(j+1)*Xkm1j{j}; ;     %x_k at the last collocation point
  
       % Add contribution to quadrature function
       J = J + B(j+1)*qj*h;
    end    
    
        % path constraints
        fIDj=F([Xkm1; Uakm1*s.ua]);%ID call
       [unitVectors, forceValues, cg] = calculateUnitVectors(fIDj, fixedPoints,Xkm1(1:2:5));
       

        g = {g{:} (forceValues-Ukm1*s.uT)/s.uT};  %equations of motion as path constraints at mesh points
        lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
        ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

    
    % Add equality constraints
    g = {g{:}, (Xk_end-Xk)};                        %equality constraint between consecutive mesh points
    lbg = [lbg; zeros(2*nq,1)];                      %lower bound for equality constraints between consecutive mesh points
    ubg = [ubg; zeros(2*nq,1)];                      %upper bound for equality constraints between consecutive mesh points
    
    %save states and controls of the current mesh point
    Ukm1=Uk;
    Uakm1=Uak;
    Xkm1=Xk;
    Xkm1j=Xkj; 
end
  
% path constraints
fIDj=F([Xkm1; Uakm1*s.ua]);        %ID call
[unitVectors, forceValues, cg] = calculateUnitVectors(fIDj, fixedPoints,Xkm1(1:2:5));
g = {g{:}, (forceValues-Ukm1*s.uT)/s.uT};  %equations of motion as path constraints at the last minus one point
lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

%States and controls at the last point
Xk = MX.sym(['X_' num2str(k+1)], nq*2);
w = {w{:}, Xk};
lbw = [lbw; final_states_lb];              %lower bound for last collocation point of x_k %%%%%%%%%%%%%%%%%%%%%%%%%
ubw = [ubw; final_states_ub];              %upper bound for last collocation point of x_k %%%%%%%%%%%%%%%%%%%%%%%%%
w0 = [w0;   final_states_ig];

Uk = MX.sym(['U_' num2str(k+1)],ncT);
w = {w{:}, Uk};  
lbw = [lbw; controls_lb];        %lower bound for u_k
ubw = [ubw; controls_ub];        %upper bound for u_k
w0 = [w0; controls_ig];          %initial guess for u_k

Uak = MX.sym(['U_a' num2str(k+1)],nq);
w = {w{:},Uak};
lbw = [lbw; accelerations_lb];     %lower bound for ua_k
ubw = [ubw; accelerations_ub];      %upper bound for ua_k
w0 = [w0; accelerations_ig];        %initial guess for ua_k

%Interpolate controls with lagrange polynomials at the collocation points of the last mesh interval
for j=1:d
    Ukj{j}=MX(ncT,1);
    Uakj{j}=MX(nq,1);
    Ukj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Ukm1 Uk]);
    Uakj{j}=LagrangePoly_CASADI(tau_root(j+1),[0 1],[Uakm1 Uak]);
end

Xk_end = D(1)*Xkm1;
for j=1:d
   % Expression for the state derivative at the collocation point
   xp = C(1,j+1)*Xkm1;
   for r=1:d
       xp = xp + C(r+1,j+1)*Xkm1j{r};
   end

   % Append collocation equations
   % defect constraints
   [fj qj]=f(Xkm1j{j},Ukj{j},Uakj{j}*s.ua);
   g = {g{:}, h*fj-xp};                 %build defect constraints
   lbg = [lbg; zeros(2*nq,1)];             %lower bounds for defect constraints
   ubg = [ubg; zeros(2*nq,1)];             %upper bounds for defect constraints

   % Add contribution to the end state
   Xk_end = Xk_end + D(j+1)*Xkm1j{j};     %x_k at the last collocation point

   % Add contribution to quadrature function
   J = J + B(j+1)*qj*h;
end  

% path constraints
fIDj=F([Xk; Uak*s.ua]);  %ID call
[unitVectors, forceValues, cg] = calculateUnitVectors(fIDj, fixedPoints,Xk(1:2:5));
g = {g{:}, (forceValues-Uk*s.uT)/s.uT};  %equations of motion as path constraints at the last minus one point
lbg = [lbg; zeros(nq,1)];             %lower bounds for path constraints
ubg = [ubg; zeros(nq,1)];             %upper bounds for path constraints

g = {g{:}, (Xk_end-Xk)};                        %equality constraint between consecutive mesh points
lbg = [lbg; zeros(2*nq,1)];                      %lower bound for equality constraints between consecutive mesh points
ubg = [ubg; zeros(2*nq,1)]; 
    
%% Create an NLP solver
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

opt.ipopt.max_iter=2000;
opt.ipopt.hessian_approximation = 'limited-memory';
solver = nlpsol('solver', 'ipopt', prob,opt);

%% Solve the NLP
sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw,...
            'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Postprocess data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get values for design variables
x_opt=[];
uT_opt=[];
ua_opt=[];
x_opt(1,1:2*nq)=w_opt(1:2*nq);
x_opt_ext=x_opt;
x_opt_ext=[x_opt_ext; reshape(w_opt((2*nq+1):(2*nq+d*2*nq)),2*nq,d)'];
uT_opt(1,1:ncT)=w_opt((2*nq+d*2*nq+1):(2*nq+d*2*nq+ncT));

ua_opt(1,1:nq)=w_opt((2*nq+d*2*nq+2+1):(2*nq+d*2*nq+2+nq));
nvarxint=2*nq+nq+ncT+2*nq*d;
for i=1:N-1
    x_opt=[x_opt; w_opt((2*nq+d*2*nq+ncT+nq+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+(i-1)*nvarxint)))'];
    uT_opt=[uT_opt; w_opt((2*nq+d*2*nq+ncT+nq+2*nq+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+ncT +(i-1)*nvarxint )))'];
    ua_opt=[ua_opt; w_opt((2*nq+d*2*nq+ncT+nq+2*nq+ncT+1 +(i-1)*nvarxint):((2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq +(i-1)*nvarxint)))'];
    x_opt_ext=[x_opt_ext; x_opt(end,:); reshape(w_opt((2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq+1 +(i-1)*nvarxint):(2*nq+d*2*nq+ncT+nq+2*nq+ncT+nq+2*nq*d +(i-1)*nvarxint)),2*nq,d)'];
end
x_opt=[x_opt; w_opt((   2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+1):((         2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq)))'];
uT_opt=[uT_opt; w_opt(( 2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+1):(     2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT))'];
ua_opt=[ua_opt; w_opt(( 2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT+1):(   2*nq+d*2*nq+ncT+nq+(N-1)*nvarxint+2*nq+ncT+nq))'];
x_opt_ext=[x_opt_ext; x_opt(end,:)];


uT_opt_ext(1:4:(N*(d+1)+1),:)=uT_opt;
ua_opt_ext(1:4:(N*(d+1)+1),:)=ua_opt;
for i=1:N
    for j=1:d
        for qi=1:nq;
            uT_opt_ext((i-1)*(d+1)+j+1,qi)=LagrangePoly(tau_root(j+1),[0 1],[uT_opt(i,qi) uT_opt(i+1,qi)]);
            ua_opt_ext((i-1)*(d+1)+j+1,qi)=LagrangePoly(tau_root(j+1),[0 1],[ua_opt(i,qi) ua_opt(i+1,qi)]);
        end
    end
end

uT_opt_ext=uT_opt_ext*s.uT;
ua_opt_ext=ua_opt_ext*s.ua;

% Gráfica de \theta (solo \theta_1, \theta_2, \theta_3)
subplot(1,2,1);
title('\theta [º]');
plot(tgrid_ext, x_opt_ext(:,1:2:3*2)*180/pi, 'LineWidth', 2);
xlim([0 1]);
xlabel('time [s]');
ylabel('\theta [º]');
legend({'\theta_1', '\theta_2', '\theta_3'}, 'Location', 'best');
hold all;
% Gráfica de \dot{\theta} (solo \dot{\theta}_1, \dot{\theta}_2, \dot{\theta}_3)
subplot(1,2,2);
title('\theta dot [rad/s]');
plot(tgrid_ext, x_opt_ext(:,2:2:3*2), 'LineWidth', 2);
xlim([0 1]);
xlabel('time [s]');
ylabel('\theta dot [rad/s]');
legend({'\dot{\theta}_1', '\dot{\theta}_2', '\dot{\theta}_3'}, 'Location', 'best');
hold all;
% Gráfica de Momentos
figure;
title('Force values');
plot(tgrid_ext, uT_opt_ext, 'LineWidth', 2);
xlim([0 1]);
xlabel('time [s]');
ylabel('M [Nm]');
legend(arrayfun(@(i) ['M_', num2str(i)], 1:size(uT_opt_ext, 2), 'UniformOutput', false), 'Location', 'best');
hold all;

% Gráfica de Valores de Fuerza (forcevalues)
% figure;
% title('Force Values');
% plot(tgrid_ext, forceValues, 'LineWidth', 2);
% xlim([0 1]);
% xlabel('time [s]');
% ylabel('Force [N]');
% legend('Force', 'Location', 'best');
% grid on;


keyboard;
%Recalculate path constraints
for i=1:size(x_opt,1)
    fIDj=F([x_opt(i,:)'; ua_opt(i,:)']);  %ID call
    fIDj_all(i,:)=full(fIDj);
    [unitVectors, forceValues, cg] = calculateUnitVectors_val(fIDj_all(i,:), fixedPoints, x_opt_ext(i,1:2:5));
    forceValues_all(i,:)=forceValues;
end

keyboard;
% Create mot file
q_str.data=[tgrid_ext' x_opt_ext(:,1:2:3*2)*180/pi x_opt_ext(:,7:2:nq*2)];
q_str.labels={'time','r1_z','r2_x','r3_y','t4_x','t5_y','t6_z'};
write_motionFile(q_str,'motion3.mot');

keyboard;

%% Evaluate sparsity of jacobian and hessian
g = vertcat(g{:});
w = vertcat(w{:});
Jac = jacobian(g, w);
figure(4)
spy(sparse(DM.ones(Jac.sparsity())))

