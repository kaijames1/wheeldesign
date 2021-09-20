% SHAPE OPTIMIZATION MAIN PROGRAM
%  for rigid-body, left-to-right rolling "wheels" over various
%   terrain types.
%  Objective function type: WORK MINIMIZATION,
%   i.e. f = sum(F.dS)
%  subject to a perimeter constraint, P_low <= P <= P_upp
%
% PATRICK KELLEY | KAI JAMES - UIUC - 2021
%
clear; clc; close all;

%**************************************************************
%** USER INPUTS ***********************************************
%**************************************************************
%**
%** ns : total number of spokes in the shape
ns = 48;
    %-- Dependencies: DO NOT ALTER! ---------------------------
    %----nv : total number of design variables; the spokes plus
    %          all the angles between them.
    %--- dQ : "default" sector angle: the interior angle in a
    %          polygon with nSpokes sides.
    nv = 2*ns;
    dQ = 2*pi/ns;
    %----------------------------------------------------------
    
  %** ANGULAR CONTROL *****************************************
  %* wedgeang : the angle a spoke can be within.              *
  %              0 < wedgeang <= pi/4, preferably. Can be a   *
  %              multiple of dQ!                              *
  %************************************************************
    wedgeang = 0.5*dQ;
    
  %** RADIAL CONTROL ******************************************
  %**** R_low : lower bound on the spoke length, [m]          *
  %**** R_upp : upper bound on the spoke length, [m]          *
  %************************************************************
    R_low = 0.1;
    R_upp = 1.2;

  %** PERIMETER CONSTRAINT CONTROL ****************************
  %  Perimeter lies between P_low*arc_c <= P <= P_upp*arc_c,  *
  %   where arc_c is the arc length over half of the terrain's*
  %   wavelength (i.e. arc length of t(0:L/2)).               *
  %**** P_low : multiplier for lower bound on the perimeter.  *
  %**** P_upp : multiplier for upper bound on the perimeter.  *
  %************************************************************
    P_low = 1.5;
    P_upp = 2.5;

  %** TERRAIN PARAMETERIZATION ********************************
  %******** A : amplitude of the terrain, [m]. Give as 0 for  *
  %              flat ground.                                 *
  %******** L : wavelength of the terrain, [m]                *
  %************************************************************
    A_T = 0;
    L_T = 3;

  %** SIMULATION PARAMETERS ***********************************
  %****** Fax : force on the axle in the vertical direction,  *
  %              [N]. Give as a positive scalar.              *
  %************************************************************
    Fax = 1e1;  % Set to 10 to keep objective function value>1

  %** OPTIMIZATION PARAMETERS *********************************
  % maxoutit : maximum number of iterations for the optimizer.*
  %************************************************************
    maxoutit = 5000;
    

%--------------------------------------------------------------
%-- COMPUTE INITIAL DEPENDENCIES
%--------------------------------------------------------------

%-- Create the matrices/vectors relating the n design variables
%    to the n+1 physical variables using a relation of the form
%                       rho = Ww*x + bw
%    where Ww is an (n+1 x n) matrix and bw is an (n+1 x 1) 
%    vector. For sensitivities relating the design variables
%    "rho" to the physical variables "x", W' suffices.
[Wb,bb] = wheelfiltmat(R_low,R_upp,ns,wedgeang);

%-- Create the matrices/vectors for converting the n+1 design
%    variables "rho" to their cartesian counterparts. bc
%    is an (nSpokes x 1) vector with each entry being the
%    cumulative angular location of the jth spoke with regards
%    to the first one at -pi/2
bc = (0:(ns-1))'*2*pi/ns - pi/2;

%-- Summation matrix for axle/pivot expansion
Summer = zeros(2,2*ns);
Summer(1,1:ns) = 1; Summer(2,ns+1:end) = 1;

%-- Generate the terrain and the arc length.
[TX,dTXdX,arc_c] = genTX(A_T,L_T);

%-- Compute the perimeter constraints.
MinP = P_low*arc_c;
MaxP = P_upp*arc_c;
    
    
%--------------------------------------------------------------
%-- RANDOM SHAPE GENERATION & FIRST COMPUTATIONS
%--------------------------------------------------------------
% Randomly initialize design variables
x0 = rand(nv,1);
tic;

% Convert to physical variables
rho0 = Wb*x0 + bb;

% Compute the work done by the initial shape
[W,dWdx] = rollsf(x0,ns,nv,Summer,Wb,bb,bc,Fax,TX,dTXdX,A_T,L_T,0);

% Compute the initial perimeter
[P,dPdrho] = chperim(rho0,ns);
% SENSITIVITY: Perimeter sensitivities are known analytically if the phys.
%  variables sent in were not convex-hulled
    dPdx = Wb'*dPdrho;
  
% Plot the shape
T = toc;
figure(2);
plotwh(rho0,ns,bc);
plotopt(arc_c,A_T,L_T);
animation(1) = getframe(); hold off;

%makepng(animation(1),'SOcircstart');

%%
%--------------------------------------------------------------
%-- INITIALIZE PARAMETERS USED IN MMA PROBLEM -----------------
%--------------------------------------------------------------    
xval = x0;                      % Initial variable values
m = 2;                          % No. of constraints
xold1   = xval;                 % Previous var. are the initial
xold2   = xval;                 % "                           "
xmin    = 0*ones(nv, 1);         % Minimum allowable values
xmax    = ones(nv, 1);      % Maximum "              "
low     = xmin;
upp     = xmax;
c_mma   = 10000*ones(m, 1);
d_mma   = 1*ones(m, 1);
a0_mma  = 1;
a_mma   = 0*ones(m, 1);
outeriter = 0;                  % Initial iteration number
kkttol  = 1.0e-2;               % KKT condition tolerance
kktnorm = kkttol + 1;           % Initial KKT gradient norm
f0val = W;                      % Initial work of the shape
df0dx = dWdx;                   % Initial gradient
fval = [MinP - P, P - MaxP];    % Constraint value
dfdx = [-dPdx, dPdx];           % Constraint gradient

% Store the previous perimeter sensitivities for oscillation damping
dfdxold1 = dfdx;


%-- Continuation parameter that limits the design search space of the
%    optimizer, as suggested in Guest 2011 for Topology Optimization
%    Heaviside Filtering. It is useful here because it ensures
%    smoothness of the final shape by severly limiting the move space
%    as this parameter undergoes continuation up to 5000
HBeta = 5;

%-- Print the table header for the iterations
fprintf('Iter.  KKT norm      f(x)      P/Pmin    %ct     %ct   HBeta\n',916,8747);
fprintf('----- ---------- ------------ -------- ------ ------ -----\n');
fprintf('%5d %10.4f %12.6f %8.4f %6.3f %6.2f %5d\n',...
    outeriter,kktnorm,W,P/arc_c,T,T/60,HBeta);

%-- History storage variable initialization
x_hist = zeros(nv,maxoutit+1);
W_hist = zeros(1,maxoutit); P_hist = zeros(1,maxoutit);
T_hist = zeros(1,maxoutit); T_hist(1) = T;
K_hist = zeros(1,maxoutit);

%--------------------------------------------------------------
%-- SOLVING THE OPTIMIZATION PROBLEM! -------------------------
%-------------------------------------------------------------- 
tic;
%%
while outeriter < maxoutit && kktnorm > kkttol
    % Increase the iteration number
        outeriter = outeriter + 1;    
        
    % Continuation of HBeta (search range for the optimizer)
        HBeta = continuation(HBeta,outeriter,HBeta,100,10e3);
        
    % The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(m,nv,outeriter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma,HBeta,dfdxold1');
    
    % Some vectors are updated:
        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;

    % Compute the new spoke and angle values
        rho = Wb*xval + bb;
        
    % Evaluate the objective function, constraints
        [W,dWdx] = rollsf(xval,ns,nv,Summer,Wb,bb,bc,Fax,TX,dTXdX,A_T,L_T,0);
        [P,dPdrho] = chperim(rho,ns);
        dPdx = Wb'*dPdrho;
        
        % To find where constraint deriv changes sign
        dfdxold1 = dfdx;

        
        
    % Update the values
        f0val = W;                      % Update the objective function value
        df0dx = dWdx;                   % Update the variable gradient      
        fval = [MinP - P, P - MaxP];    % Check the perimeter constraint value
        dfdx = [-dPdx, dPdx];           % Update perimeter constraint gradient
        
    % The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,nv,xmma,ymma,zmma,lam,xsi,eta,mu,zet,...
        s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);
    
    % Output Information
        T = T + toc;                    % Keep track of the total time
        fprintf('%5d %10.4f %12.6f %8.4f %6.3f %6.2f %5d\n',outeriter,kktnorm,W,P/arc_c,toc,T/60,HBeta);
        tic;                            % Reset the timer
        
    % Plot Wheel Design
        figure(2);
        plotwh(rho,ns,bc);
        plotopt(MinP/2,A_T,L_T);
        title(['Iteration: ',num2str(outeriter)]); axis off;
        getframe(); hold off;
%         animation(outeriter+1) = getframe(2); hold off;
        
    % Save the current status to the history arrays    
        x_hist(:, outeriter) = xval;
        W_hist(outeriter) = W;
        P_hist(outeriter) = P;
        T_hist(outeriter) = T;
        K_hist(outeriter) = kktnorm;
    % Repeat...
end

%% ------------------------------------------------------------
%-- POST-PROCESSING -------------------------------------------
%--------------------------------------------------------------

%-- Inform the user of the optimizer's reason for stopping.
if kktnorm < kkttol
    fprintf('\nKKT conditions satisfied to within tolerance!\n');
else
    fprintf('\nMaximum optimizer iterations reached.\n');
end

%% Plot objective function history
    figure(3); cla
    semilogy(W_hist,'k');
    set(gca,'FontName','Times New Roman');
    xlim([0 enditer]);
    title('Objective function value across iterations');
    xlabel('Iteration');
    ylabel('Work, f(\bf{x}\rm)');
    grid on;
    
%% Plot iteration time history    
    figure(4); cla
    plot(diff(T_hist),'k');
    set(gca,'FontName','Times New Roman');
    xlim([0 outeriter-2]);
    title('Iteration time across iterations');
    xlabel('Iteration');
    ylabel('Time/iteration, [s]');
    grid on;
%% Plot KKTnorm time history    
    figure(5); cla
    semilogy(P_hist,'k');
    set(gca,'FontName','Times New Roman');
    xlim([0 outeriter-2]);
    title('KKT norm across iterations');
    xlabel('Iteration');
    ylabel('KKT norm');
    grid on;
%%    
% Show the optimized shape rolling over the terrain (and get its animation)
    figure(10);
    [~,~,animroll] = rollsf(xval,ns,nv,Summer,Wb,bb,bc,Fax,TX,dTXdX,A_T,L_T,1);

    
%% SAVING RESULTS    
% makegif(animation(:),'newcirc');
% 
% if A_T == 0
%     tystr = 'circle';
% else
%     tystr = 'ellipse';
% end
% 
% enditer = outeriter;
% 
% for i = [1:100:round(enditer-100,-2)+1,enditer]
%     
%     rho = Wb*x_hist(:,i) + bb;
%     figure(8);
%     plotwh(rho,ns,bc);
%     plotopt(MinP/2,A_T,L_T);
%     
%     makepng(getframe(8),['SO_',num2str(tystr),'_iter',num2str(i-1)]);
%     hold off;
%     
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [TX,dTXdX,arc_c] = genTX(A_T,L_T)
%GENTX     Generate the terrain profile, its derivative WRT X, and
%           compute the arc length for the perimter constraint
%
%  INPUTS:
%   A_T.....Amplitude of the ground, [m]
%   L_T.....Wavelength of the terrain, [m]
%  OUTPUTS:
%   TX......Terrain function T(X)
%   dTXdX...Slope of the terrain at each point X
%   arc_c...Arc length of the terrain over half a wavelength, [m]
%
% PATRICK KELLEY - UIUC - 2021

    %** The FLAT ground case is super easy, barely an inconvenience
    if A_T == 0
        TX = @(X) 0;
        dTXdX = @(X) 0;
        arc_c = L_T/2;
    %** The SINUSOIDAL terrain case
    else
        TX = @(X) A_T.*sin(pi.*X./L_T).^2;
        dTXdX = @(X) 2*A_T*sin(pi*X/L_T)*cos(pi*X/L_T)*pi/L_T;
        dx = 1e-3;
        Gx = 0:dx:L_T/2;
        arc_c = trapz(sqrt(1+(A_T*pi/L_T*sin(2*pi/L_T*Gx)).^2))*dx;
    end
    %** Perhaps more functions in the future!
end


function [W,b] = wheelfiltmat(R_low,R_upp,ns,wedge)
%WHEELFILTMAT   Create the matrix for converting the design
%                variables to the physical variables in a 
%                wheel shape optimization problem.
%
%   [W,b] = WHEELFILTMAT(R_low,R_upp,nSpokes,wedge) creates
%    the necessary matrices and vectors so the conversion of
%    wheel design variables "x" to their physical realizations
%    "rho" is of the form
%                           rho = W*x + b
%    which neatly give the sensitivities of rho WRT x as W'.
%   The form and values in W and b depend on the ratio of the
%    wedge over the default angle dQ = 2*pi/nSpokes. If wedge
%    < dQ, then a simple linear relationship is given. However,
%    if wedge >= dQ, a more complicated form is supplied as to
%    ensure spokes cannot overlap or get too close.
%
%  INPUTS:
%   R_low.....Lower bound on the spoke length, [m]
%   R_upp.....Upper bound on the spoke length, [m]
%   ns........Number of spokes in the wheel.
%   wedge.....Sector angle [rad] within which a spoke can be
%              present
%  OUTPUTS:
%   W.........2*nSpokes x (2*nSpokes - 1) matrix
%   b.........2*nSpokes x 1 constant vector
%
% PATRICK KELLEY - UIUC - 2020
    
    % CREATE THE SUB-MATRIX AND CONSTANT VECTOR FOR RADII
    Wr = eye(ns)*(R_upp - R_low);
    br = ones(ns,1)*(R_low);
     
    % CREATE THE SUB-MATRIX AND CONSTANT VECTOR FOR ANGLES
    dQ = 2*pi/ns;      % Default sector angle
    
    if wedge < dQ
        % If the wedge angle is smaller than the default sector,
        %  then the angle varies linearly between the lower bound
        %  and upper bound, with no dependence on adjacent spokes
        %  aside from creating the differences as Q requires
        Wq = eye(ns)*wedge;     

        % Constant vector
        bq = ones(ns,1)*(-wedge/2);
    else
        % To be included, if necessary
    end
    
    % CREATE THE GLOBAL ARRAYS
    W = blkdiag(Wr,Wq);
    b = [br;bq];
end


function [Xi,a,dXidrho,dadrho] = physvar2Xi(rho,ns,bc)
%PHYSVAR2Xi   Convert radial/angular physical variables
%             into complex numbers using polar form.
%
% INPUTS:
%   rho.......(2ns x 1) physical variables. First ns are the
%              spoke lengths, next ns are converted angles [rad]
%   ns........number of spokes in the shape
%   bc........(ns x 1) vector for shifting the angles
%
% OUTPUTS:
%   Xi........(2ns x 1) locations of the spoke tips, [X;Y]
%   a.........(2 x 1) axle location
%   dXidrho...(2ns x 2ns) sensitivity matrix of all coordinates
%   dadrho....(2 x 2ns) sensitivity matrix for the axle location
%
% PATRICK KELLEY - UIUC 2021

    % Prepare
    Xi = zeros(2*ns,1);
    
    % Extract
    R = rho(1:ns);
    Q = rho(ns+1:end);
    
    % Axle
    a = [R(1)*cos(Q(1) + bc(1));R(1)*sin(Q(1) + bc(1))];
    
    % All tips
    Xi(1:ns)     = R.*cos(Q + bc) - R(1)*cos(Q(1) + bc(1));
    Xi(ns+1:end) = R.*sin(Q + bc) - R(1)*sin(Q(1) + bc(1));
    
    % SENSITIVITY of tip coordinates WRT rho
    if nargout >= 3
        % Prepare
        dXidrho = zeros(2*ns,2*ns);
        % Radial sensitivities
        dXidrho(1:ns,1:ns) = diag(cos(Q + bc));
        dXidrho(ns+1:end,1:ns) = diag(sin(Q + bc));
        dXidrho(1:ns,1) = dXidrho(1:ns,1) - cos(Q(1) + bc(1));
        dXidrho(ns+1:end,1) = dXidrho(ns+1:end,1) - sin(Q(1) + bc(1));
        % Angular sensitiviites
        dXidrho(1:ns,ns+1:end) = diag(-R.*sin(Q + bc));
        dXidrho(ns+1:end,ns+1:end) = diag(R.*cos(Q + bc));
        dXidrho(1:ns,ns+1) = dXidrho(1:ns,ns+1) + R(1)*sin(Q(1) + bc(1));
        dXidrho(ns+1:end,ns+1) = dXidrho(ns+1:end,ns+1) - R(1)*cos(Q(1) + bc(1));
    end
    
    % SENSITIVITY of axle coordinates WRT rho
    if nargout == 4
        % Prepare
        dadrho = zeros(2,2*ns);
        % X-location sensitivity
        dadrho(1,1) = cos(Q(1) + bc(1));
        dadrho(1,ns+1) = -R(1)*sin(Q(1) + bc(1));
        % Y-location sensitivity
        dadrho(2,1) = sin(Q(1) + bc(1));
        dadrho(2,ns+1) = R(1)*cos(Q(1) + bc(1));
    end
    
end


function [Wtot,dWtotdx,anim] = rollsf(x,ns,nv,Summer,Wb,bb,bc,F,TX,dTdX,A_T,L_T,shoplot)
%ROLLSF  Compute the work done by a rolling shape over sinusoidal terrain.
%   
%   [Wtot,dWtotdx] = ROLLSF(x,ns,nv,Summer,Wb,bb,Wc,bc,F,L_T,shoplot)
%
%  INPUTS:
%   x..........design variable vector, (nv x 1)
%   ns ... ... number of spokes
%   nv.........number of design variables (2ns)
%   Summer ... (2 x 2ns) matrix for summing X-coords with
%               the first row, Y-coords with the second
%   Wb.........design 2 physical matrix, (2ns x nv)
%   bb ... ... design 2 physical offset vector (2ns x 1)
%   bc ... ... cartesian offset vector
%   F..........Force acting downward through the axle, [N]
%   TX.........terrain function of X-location
%   dTdX.. ... derivative function of the the terrain wrt X-location
%   A_T........amplitude of the ground, [m]
%   L_T... ... Wavelength of the ground, [m]
%   shoplot....Plotting toggle: 0 if no plots, 1 if yes
%
%  OUTPUTS:
%   Wtot.......Total work done by the shape's axle over a roll
%   dWtotdx... (nv x 1) gradient vector of the work function wrt x
%
% PATRICK KELLEY - KAI JAMES - UIUC 2021

% NOTE: [FDV] tag means sensitivity was verified by Finite Difference
%       [ANA] means sens is an analytic relation verified by inspection
%========================================================================

    % Total number of rotation steps: have each spoke tip be the pivot
    %  two times so that the first spoke gets sensitivities no just from
    %  its initial position
    nr = 2*ns;
    
    % Convert from DESIGN to PHYSICAL variables
    rho = Wb*x + bb;    % (2ns x 1)
    % SENSITIVITY: physical variables wrt x [ANA]
        drhodx = Wb';   % (nv x 2ns)
    % Convert to cartesian coordinates; get sensitivity
    [Xi0,Xax,dXidrho,dXaxdrho] = physvar2Xi(rho,ns,bc);
    % SENSITIVITY: Cartesian coords WRT design vars (2ns x 2nv) [FDV]
        dXidx = zeros(2*ns,nv,nr+1);
        dXidx(:,:,1) = dXidrho*drhodx';

    %-- INITIALIZATIONS --------------------------------------
        W = zeros(1,nr+1);      % Work at each step in 'time'
        dWdx = zeros(nr+1,nv);  % sensitivity of work over 'time'
        Xi = zeros(2*ns,nr+1);  % Cartesian locations of spoke tips over 'time'
        Xi(:,1) = Xi0;          % Initial cartesian coords
        a = zeros(2,nr+1);      % axle locations over 'time'
        p = zeros(2,nr+1);      % pivot coords over 'time'
        a(:,1) = -Xax;               % initial axle location
        % SENSITIVITY: Axle WRT x (2 x 2ns) [ANA]
            dadx = zeros(2,nv,nr+1);
            dadx(:,:,1) = -dXaxdrho*drhodx;
        % SENSITIVITY: Pivot WRT x (nv x 2)
            dpdx = zeros(2,nv,nr+1);
    %-- END OF INITIALIZATIONS -------------------------------
    
    %++ Show animation of shape rolling over the terrain, if desired ++++++
    if shoplot == 1
        plot(0:1e-2:3*L_T,TX(0:1e-2:3*L_T),'k'); hold on; 
        grid on; plot(p(1,1),p(2,1),'r*'); plotwheel(ns,Xi(:,1),a(:,1));
        xlim([-L_T/2,3*L_T]); ylim([0,L_T]); title(['j = ',num2str(0)]);
        anim(1) = getframe(); hold off;
    else
        anim = 0;
    end
    
    %### BEGIN THE SIMULATION LOOP ########################################
    k = 1; 
    for j = 1:nr
        
        % Index of next spoke
        if k < ns; k = k+1;
        else; k = 1; end
        
        % Adjusted index of current spoke, for when j > ns
        if j <= ns; m = j;
        else;  m = j - ns; end
        
        % Length of the segment between jth (pivot) spoke and the next one
        [rj,drjdrho] = perimseg(rho,m,k,ns);
        % SENSITIVITY: segment WRT rho, then x [FDV]
            drjdx = (drhodx*drjdrho)';
        
        % Next pivot point
        [p(:,j+1),dpdx(:,:,j+1)] = nextpivot(A_T,nv,p(:,j),dpdx(:,:,j),rj,drjdx,TX,dTdX);
             
        %++ OPTIONAL PLOTTING (of next pivot location) ++
        if shoplot == 1; hold on; plot(p(1,j+1),p(2,j+1),'bo'); hold off; end
        
        %-- ANGLE CALCULATION ---------------------------------------------
        %** Calculate angle to rotate shape
        % 'beta' is the angle, from y=0, between the current pivot location
        %   and the next one
        [beta,dbeta] = arctan2(p(2,j+1) - p(2,j),p(1,j+1) - p(1,j));
        % SENSITIVITY: angle between pivots
            dbetadx = dbeta'*[(dpdx(1,:,j+1) - dpdx(1,:,j));(dpdx(2,:,j+1)-dpdx(2,:,j))];
        % 'gamm' is the angle between the current pivot loc. and the
        %   current location of the next spoke tip
        [gamm,dgamm] = arctan2(Xi(ns+k,j) - p(2,j), Xi(k,j) - p(1,j));
        % SENSITIVITY: angle between pivot and next spoke
            dgammdx = dgamm'*[(dXidx(k,:,j) - dpdx(1,:,j));(dXidx(ns+k,:,j)-dpdx(2,:,j))];
        % Rotation angle 'alph' is the difference between gamm (pivot to
        %   spoke tip) and beta (pivot to next pivot), and thus the angle
        %   the wheel needs to rotate
        alph = gamm - beta;
        % SENSITIVITY:  rotation angle wrt x
            dalphdx = dgammdx - dbetadx;
        
        % Rotation matrix and its sensitivities [FDV]
        [Rbar,dRbardQ] = Rbarmat(alph,ns);
        
        %-- UPDATE THE SHAPE'S LOCATION -----------------------------------
        % Updated cartesian coordinates of the shape
        Xi(:,j+1) = Rbar*Xi(:,j) + (eye(2*ns) - Rbar)*Summer'*p(:,j);
        % SENSITIVITY: Updated shape location wrt x [FVD]
            dXidx(:,:,j+1) = (dRbardQ*Xi(:,j))*dalphdx ...
                + Rbar*dXidx(:,:,j) ...
                + Summer'*dpdx(:,:,j) ...
                - dRbardQ*Summer'*p(:,j)*dalphdx ...
                - Rbar*Summer'*dpdx(:,:,j);
        % Updated axle location
        a(:,j+1) = [Rbar(1,1),Rbar(1,ns+1);Rbar(ns+1,1),Rbar(ns+1,ns+1)]*...
                    (a(:,j) - p(:,j)) + p(:,j);
        % SENSITIVITY: axle wrt new shape coords [FVD]
              dadx(:,:,j+1) = [Rbar(1,1),Rbar(1,ns+1);Rbar(ns+1,1),Rbar(ns+1,ns+1)]*...
                    (dadx(:,:,j) - dpdx(:,:,j)) + ...
                    [dRbardQ(1,1),dRbardQ(1,ns+1);dRbardQ(ns+1,1),dRbardQ(ns+1,ns+1)]*...
                   (a(:,j) - p(:,j))*dalphdx + dpdx(:,:,j);
        
        %-- WORK CALCULATION ----------------------------------------------
        % 1st cond: Moved ahead but the pivot spoke didn't pass the y-axis
        % 2nd cond: Moved backwards but didn't cross back through y-axis
            if (a(1,j+1) < p(1,j) && a(1,j+1)>a(1,j)) ...
                    || ((a(1,j+1) < a(1,j)) && (a(1,j+1) >= p(1,j)))
                W(j) = (a(2,j+1) - a(2,j));
                % SENSITIVITY: work wrt axle location wrt x [ANA]
                    dWdx(j,:) = (dadx(2,:,j+1) - dadx(2,:,j));                
            % Otherwise, move through the y-axis one way or the other, and
            %  therefrore a portion of the axle height is irrelevant
            else
                W(j) = (p(2,j) + rho(m) - a(2,j));
                % SENSITIVITY: work wrt axle location wrt x [ANA]
                    dWdx(j,:) = (dpdx(2,:,j) - dadx(2,:,j));
                    dWdx(j,m) = dWdx(j,m) + drhodx(m,m);
            end
        %-- END OF WORK CALCULATION ---------------------------------------
        
        %++ Show animation of shape rolling over the terrain, if desired ++
        if shoplot == 1
            plot(0:1e-2:3*L_T,TX(0:1e-2:3*L_T),'k'); hold on; 
            grid on; plot(p(1,j+1),p(2,j+1),'r*'); plotwheel(ns,Xi(:,j+1),a(:,j+1));
            plot(a(1,1:j+1),a(2,1:j+1),'g');
            xlim([-L_T/2,3*L_T]); ylim([0,L_T]); 
            title(['j = ',num2str(j)]); anim(j+1) = getframe(); hold off;
            %pause(0.15);
        end
    end %### End of simulation loop #######################################
    
    % Sum everything up
    Wtot = F*sum(W); 
    dWtotdx = F*sum(dWdx)';
end


function [rj,drjdrho] = perimseg(rho,j,k,ns)
    dQ = 2*pi/ns;
    %Perimeter segment length
    rj = sqrt(rho(j)^2 + rho(k)^2 - 2*rho(j)*rho(k)*cos(dQ - rho(ns+j) + rho(ns+k)));
    % SENSITIVITY
    if nargout > 1
        drjdrho = zeros(size(rho));
        drjdrho(j) = 1/rj*(rho(j) - rho(k)*cos(dQ - rho(ns+j) + rho(ns+k)));
        drjdrho(k) = 1/rj*(rho(k) - rho(j)*cos(dQ - rho(ns+j) + rho(ns+k)));
        drjdrho(ns+j) = -1/rj*(rho(j)*rho(k)*sin(dQ - rho(ns+j) + rho(ns+k)));
        drjdrho(ns+k) = 1/rj*(rho(j)*rho(k)*sin(dQ - rho(ns+j) + rho(ns+k)));
    end
end


function [p,dpdx1] = nextpivot(A_T,nv,pj,dpdxj,rj,drjdx,T,dTdX)

    % Flat ground case
    if A_T == 0
        p = pj + [rj;0];
        % SENSITIVITY: new pivot point wrt x [ANA]
            dpdx1 = dpdxj + [drjdx;zeros(1,nv)];
    else
        tol = 1e-3;
        pik = pj(1) + rj;
        % SENSITIVITY: initial pivot guess WRT x [ANA]
            dpikdx = [1 0; 0 0]*dpdxj + drjdx;
        % Big initial change to get the loop going (will be overwritten)
        Dpi = 100;

        % Calculate the x-location of the next pivot point
        while abs(Dpi) > tol
            fhh = T(pik)^2 + pik^2 - 2*(pj(1)*pik+pj(2)*T(pik)) + pj(1)^2 + pj(2)^2 - rj^2;
            fhhprime = 2*(T(pik)*dTdX(pik) + pik - (pj(1) + pj(2)*dTdX(pik)));
            % SENSITIVITY: function wrt x
                dfhhdx = fhhprime*dpikdx + 2*((pj(1)-pik)*dpdxj(1,:) + (pj(2)-T(pik))*dpdxj(2,:) - rj*drjdx);
            % SENSITIVITY: slope function wrt x
                dfhhprimedx = 2*(((T(pik)-pj(2))*dTdX(pik)+1)*dpikdx - dpdxj(1,:) - dTdX(pik)*dpdxj(2,:));
            Dpi = -fhh/fhhprime;
            pik = pik + Dpi;
            % SENSITIVITY: pivot point wrt x
                dpikdx = dpikdx - dfhhdx/fhhprime + fhh*dfhhprimedx/fhhprime^2;
        end
        
        % Next pivot point
        p(1) = pik;
        p(2) = T(pik);
        % SENSITIVITY: new pivot point wrt x
            dpdx1 = [1 0;0 dTdX(pik)]*dpikdx;
    end
end


function [t, dtdxy] = arctan2(y,x)
%ARCTAN2    The atan2 function, but with an analytic
%            gradient as an optional second output
%
%   [t,dtdxy] = ARCTAN2(y,x)
%
%  INPUTS:
%   y.......Y-coordinate
%   x.......X-coordinate
%  OUTPUTS:
%   t.......signed angle between -pi and pi; [rad]
%   dtdxy...gradient vector, [dtdx;dtdy]
%
% PATRICK KELLEY - 2021

    t = atan2(y,x);

    if nargout > 1
        nxy = x.^2 + y.^2;
        dtdxy = [-y./nxy; x./nxy];
    end
end


function [Rbar,dRbardQ] = Rbarmat(Q,ns)
    % Create rotation matrix blocks
    Cblock = eye(ns)*cos(Q);
    Sblock = eye(ns)*sin(Q);
    % Assemble
    Rbar = [Cblock,Sblock;-Sblock,Cblock];
    
    % SENSITIVITY of rotation matrix wrt the rotation angle
    if nargout > 1
        dCblock = -eye(ns)*sin(Q);
        dSblock = eye(ns)*cos(Q);
        dRbardQ = [dCblock,dSblock;-dSblock,dCblock];
    end
end


function plotwheel(ns,Xi,a)
    plot([Xi(1:ns); Xi(1)],[Xi(ns+1:end); Xi(ns+1)],'-ok');
    axis equal; hold on; grid on;
    text(a(1),a(2),'@','HorizontalAlignment','center','VerticalAlignment','middle');
    plot([a(1),Xi(1)],[a(2),Xi(ns+1)],'r');
end


function [P,dPdC] = chperim(C,ns)
%CHPERIM   Compute the perimeter of a convex hull defined by
%           radial (R) and angular (Q) variables C = [R;Q].
%
%  INPUTS:
%   C.........design variable vector, R [m] and Q [rad], with
%              R having passed through a smooth convex hull
%   ns........number of spokes in the wheel
%
%  OUTPUTS:
%   P.........Perimeter length, [m]
%   (dPdC)....Optional: the sensitivity vector of the perimeter
%              WRT to the convex hull design variables
%
% PATRICK KELLEY - UIUC - 2020
    
    % Extract the radial and angular variables
    R = C(1:ns,1);
    Q = C(ns+1:end,1);
    dQ = 2*pi/ns;
    % Preallocate the sensitivities if needed
    if nargout > 1
        dPdC = zeros(2*ns,1);
    end
    
    % COMPUTE THE PERIMETER
    P = 0;
    
    for j = 1:ns

        if j+1 > ns
            k = 1;
        else
            k = j+1;
        end

        % Calculate the length of the current segment
        Pp = sqrt(R(j)^2 + R(k)^2 - 2*R(j)*R(k)*cos(dQ-Q(j)+Q(k)));
        % Add it to the running total
        P = P + Pp;
        
        % SENSITIVITY
        if nargout > 1
            % Radial sensitivities
            dPdC(j) = dPdC(j) + 1/Pp*(R(j) - R(k)*cos(dQ - Q(j) + Q(k)));
            dPdC(k) = dPdC(k) + 1/Pp*(R(k) - R(j)*cos(dQ - Q(j) + Q(k)));
            % Angular 
            dPdC(ns+j) = dPdC(ns+j) - 1/Pp*(R(j)*R(k)*sin(dQ - Q(j) + Q(k)));
            dPdC(ns+k) = dPdC(ns+k) + 1/Pp*(R(j)*R(k)*sin(dQ - Q(j) + Q(k)));
        end
    end
end


%----------------------------------------------------------
% MMA FUNCTIONS
%----------------------------------------------------------


%-------------------------------------------------------
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,H_beta,dfdxold1)
%
%    Version September 2007 (and a small change August 2008)
%
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics, SE-10044 Stockholm, Sweden.
%
%    This function mmasub performs one MMA-iteration, aimed at
%    solving the nonlinear programming problem:
%         
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmin_j <= x_j <= xmax_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%  iter  = Current iteration number ( =1 the first time mmasub is called).
%  xval  = Column vector with the current values of the variables x_j.
%  xmin  = Column vector with the lower bounds for the variables x_j.
%  xmax  = Column vector with the upper bounds for the variables x_j.
%  xold1 = xval, one iteration ago (provided that iter>1).
%  xold2 = xval, two iterations ago (provided that iter>2).
%  f0val = The value of the objective function f_0 at xval.
%  df0dx = Column vector with the derivatives of the objective function
%          f_0 with respect to the variables x_j, calculated at xval.
%  fval  = Column vector with the values of the constraint functions f_i,
%          calculated at xval.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at xval.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%  low   = Column vector with the lower asymptotes from the previous
%          iteration (provided that iter>1).
%  upp   = Column vector with the upper asymptotes from the previous
%          iteration (provided that iter>1).
%  a0    = The constants a_0 in the term a_0*z.
%  a     = Column vector with the constants a_i in the terms a_i*z.
%  c     = Column vector with the constants c_i in the terms c_i*y_i.
%  d     = Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%     
%*** OUTPUT:
%
%  xmma  = Column vector with the optimal values of the variables x_j
%          in the current MMA subproblem.
%  ymma  = Column vector with the optimal values of the variables y_i
%          in the current MMA subproblem.
%  zmma  = Scalar with the optimal value of the variable z
%          in the current MMA subproblem.
%  lam   = Lagrange multipliers for the m general MMA constraints.
%  xsi   = Lagrange multipliers for the n constraints alfa_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - beta_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general MMA constraints.
%  low   = Column vector with the lower asymptotes, calculated and used
%          in the current MMA subproblem.
%  upp   = Column vector with the upper asymptotes, calculated and used
%          in the current MMA subproblem.
%
%epsimin = sqrt(m+n)*10^(-9);

    if nargin < 19
        H_beta = 0;
    end

    epsimin = 10^(-7);
    raa0 = 0.00001;
    albefa = 0.6;
    asyinit = 0.5/(1+H_beta);   % Modified by PK
    asyincr = 1.1;              % "
    asydecr = 0.9;              % "
    eeen = ones(n,1);
    eeem = ones(m,1);
%     zeron = zeros(n,1);

    % Calculation of the asymptotes low and upp :
    if iter < 2.5
        low = xval - asyinit*(xmax-xmin);
        upp = xval + asyinit*(xmax-xmin);
    else
        zzz = (xval-xold1).*(xold1-xold2);
        factor = eeen;
        %   factor(find(zzz > 0)) = asyincr;
        %   factor(find(zzz < 0)) = asydecr;
        factor((zzz > 0)) = asyincr;
        factor((zzz < 0)) = asydecr;
        low = xval - factor.*(xold1 - low);
        upp = xval + factor.*(upp - xold1);
        lowmin = xval - 10*(xmax-xmin);
        %   lowmax = xval - 0.01*(xmax-xmin);
        %   uppmin = xval + 0.01*(xmax-xmin);
        lowmax = xval - min(0.01,1/H_beta)*(xmax-xmin);     % Modified by PK
        uppmin = xval + min(0.01,1/H_beta)*(xmax-xmin);     % "
        uppmax = xval + 10*(xmax-xmin);                     % ", was 10 instead of 1
        low = max(low,lowmin);
        low = min(low,lowmax);
        upp = min(upp,uppmax);
        upp = max(upp,uppmin);
    end

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PK %
    Dx = xval-xold1;
    xxoavg = (xval+xold1)./2;
    A = -10*eeen;
    B = 10*eeen;
    dfdxsignI = find((dfdx(1,:).*dfdxold1(1,:)) <= 0);
    % Find indices where Dx >, < 0
    Dxg0I = find(Dx > 0);
    Dxl0I = find(Dx < 0);
    Aind = intersect(dfdxsignI,Dxg0I);
    Bind = intersect(dfdxsignI,Dxl0I);
    A(Aind) = xxoavg(Aind);
    B(Bind) = xxoavg(Bind);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Calculation of the bounds alfa and beta :
    zzz = low + albefa*(xval-low);
    alfa = max(zzz,xmin);
    zzz = upp - albefa*(upp-xval);
    beta = min(zzz,xmax);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PK %
    alfa = max(alfa,A);
    beta = min(beta,B);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

     
    % Calculations of p0, q0, P, Q and b.

    xmami = xmax-xmin;
    xmamieps = 0.00001*eeen;
    xmami = max(xmami,xmamieps);
    xmamiinv = eeen./xmami;
    ux1 = upp-xval;
    ux2 = ux1.*ux1;
    xl1 = xval-low;
    xl2 = xl1.*xl1;
    uxinv = eeen./ux1;
    xlinv = eeen./xl1;
    %
%     p0 = zeron;
%     q0 = zeron;
    p0 = max(df0dx,0);
    q0 = max(-df0dx,0);
    %p0(find(df0dx > 0)) = df0dx(find(df0dx > 0));
    %q0(find(df0dx < 0)) = -df0dx(find(df0dx < 0));
    pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
    p0 = p0 + pq0;
    q0 = q0 + pq0;
    p0 = p0.*ux2;
    q0 = q0.*xl2;
    %
%     P = sparse(m,n);
%     Q = sparse(m,n);
    P = max(dfdx,0);
    Q = max(-dfdx,0);
    %P(find(dfdx > 0)) = dfdx(find(dfdx > 0));
    %Q(find(dfdx < 0)) = -dfdx(find(dfdx < 0));
    PQ = 0.001*(P + Q) + raa0*eeem*xmamiinv';
    P = P + PQ;
    Q = Q + PQ;
    P = P * spdiags(ux2,0,n,n);
    Q = Q * spdiags(xl2,0,n,n);
    b = P*uxinv + Q*xlinv - fval ;
    %
    %%% Solving the subproblem by a primal-dual Newton method
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s] = ...
    subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d);
end


%-------------------------------------------------------------
%    This is the file subsolv.m
%
%    Version Dec 2006.
%    Krister Svanberg <krille@math.kth.se>
%    Department of Mathematics, KTH,
%    SE-10044 Stockholm, Sweden.
%
function [xmma,ymma,zmma,lamma,xsimma,etamma,mumma,zetmma,smma] = ...
subsolv(m,n,epsimin,low,upp,alfa,beta,p0,q0,P,Q,a0,a,b,c,d)
%
% This function subsolv solves the MMA subproblem:
%         
% minimize   SUM[ p0j/(uppj-xj) + q0j/(xj-lowj) ] + a0*z +
%          + SUM[ ci*yi + 0.5*di*(yi)^2 ],
%
% subject to SUM[ pij/(uppj-xj) + qij/(xj-lowj) ] - ai*z - yi <= bi,
%            alfaj <=  xj <=  betaj,  yi >= 0,  z >= 0.
%        
% Input:  m, n, low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d.
% Output: xmma,ymma,zmma, slack variables and Lagrange multiplers.
%
    een = ones(n,1);
    eem = ones(m,1);
    epsi = 1;
    epsvecn = epsi*een;
    epsvecm = epsi*eem;
    x = 0.5*(alfa+beta);
    y = eem;
    z = 1;
    lam = eem;
    xsi = een./(x-alfa);
    xsi = max(xsi,een);
    eta = een./(beta-x);
    eta = max(eta,een);
    mu  = max(eem,0.5*c);
    zet = 1;
    s = eem;
    itera = 0;
    while epsi > epsimin
      epsvecn = epsi*een;
      epsvecm = epsi*eem;
      ux1 = upp-x;
      xl1 = x-low;
      ux2 = ux1.*ux1;
      xl2 = xl1.*xl1;
      uxinv1 = een./ux1;
      xlinv1 = een./xl1;
      plam = p0 + P'*lam ;
      qlam = q0 + Q'*lam ;
      gvec = P*uxinv1 + Q*xlinv1;
      dpsidx = plam./ux2 - qlam./xl2 ;
      rex = dpsidx - xsi + eta;
      rey = c + d.*y - mu - lam;
      rez = a0 - zet - a'*lam;
      relam = gvec - a*z - y + s - b;
      rexsi = xsi.*(x-alfa) - epsvecn;
      reeta = eta.*(beta-x) - epsvecn;
      remu = mu.*y - epsvecm;
      rezet = zet*z - epsi;
      res = lam.*s - epsvecm;
      residu1 = [rex' rey' rez]';
      residu2 = [relam' rexsi' reeta' remu' rezet res']';
      residu = [residu1' residu2']';
      residunorm = sqrt(residu'*residu);
      residumax = max(abs(residu));
      ittt = 0;
      while residumax > 0.9*epsi && ittt < 200
        ittt=ittt + 1;
        itera=itera + 1;
        ux1 = upp-x;
        xl1 = x-low;
        ux2 = ux1.*ux1;
        xl2 = xl1.*xl1;
        ux3 = ux1.*ux2;
        xl3 = xl1.*xl2;
        uxinv1 = een./ux1;
        xlinv1 = een./xl1;
        uxinv2 = een./ux2;
        xlinv2 = een./xl2;
        plam = p0 + P'*lam ;
        qlam = q0 + Q'*lam ;
        gvec = P*uxinv1 + Q*xlinv1;
        GG = P*spdiags(uxinv2,0,n,n) - Q*spdiags(xlinv2,0,n,n);
        dpsidx = plam./ux2 - qlam./xl2 ;
        delx = dpsidx - epsvecn./(x-alfa) + epsvecn./(beta-x);
        dely = c + d.*y - lam - epsvecm./y;
        delz = a0 - a'*lam - epsi/z;
        dellam = gvec - a*z - y - b + epsvecm./lam;
        diagx = plam./ux3 + qlam./xl3;
        diagx = 2*diagx + xsi./(x-alfa) + eta./(beta-x);
        diagxinv = een./diagx;
        diagy = d + mu./y;
        diagyinv = eem./diagy;
        diaglam = s./lam;
        diaglamyi = diaglam+diagyinv;
        if m < n
          blam = dellam + dely./diagy - GG*(delx./diagx);
          bb = [blam' delz]';
          Alam = spdiags(diaglamyi,0,m,m) + GG*spdiags(diagxinv,0,n,n)*GG';
          AA = [Alam     a
                a'    -zet/z ];
          solut = AA\bb;
          dlam = solut(1:m);
          dz = solut(m+1);
          dx = -delx./diagx - (GG'*dlam)./diagx;
        else
          diaglamyiinv = eem./diaglamyi;
          dellamyi = dellam + dely./diagy;
          Axx = spdiags(diagx,0,n,n) + GG'*spdiags(diaglamyiinv,0,m,m)*GG;
          azz = zet/z + a'*(a./diaglamyi);
          axz = -GG'*(a./diaglamyi);
          bx = delx + GG'*(dellamyi./diaglamyi);
          bz  = delz - a'*(dellamyi./diaglamyi);
          AA = [Axx   axz
                axz'  azz ];
          bb = [-bx' -bz]';
          solut = AA\bb;
          dx  = solut(1:n);
          dz = solut(n+1);
          dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        end
    %
        dy = -dely./diagy + dlam./diagy;
        dxsi = -xsi + epsvecn./(x-alfa) - (xsi.*dx)./(x-alfa);
        deta = -eta + epsvecn./(beta-x) + (eta.*dx)./(beta-x);
        dmu  = -mu + epsvecm./y - (mu.*dy)./y;
        dzet = -zet + epsi/z - zet*dz/z;
        ds   = -s + epsvecm./lam - (s.*dlam)./lam;
        xx  = [ y'  z  lam'  xsi'  eta'  mu'  zet  s']';
        dxx = [dy' dz dlam' dxsi' deta' dmu' dzet ds']';
    %    
        stepxx = -1.01*dxx./xx;
        stmxx  = max(stepxx);
        stepalfa = -1.01*dx./(x-alfa);
        stmalfa = max(stepalfa);
        stepbeta = 1.01*dx./(beta-x);
        stmbeta = max(stepbeta);
        stmalbe  = max(stmalfa,stmbeta);
        stmalbexx = max(stmalbe,stmxx);
        stminv = max(stmalbexx,1);
        steg = 1/stminv;
    %
        xold   =   x;
        yold   =   y;
        zold   =   z;
        lamold =  lam;
        xsiold =  xsi;
        etaold =  eta;
        muold  =  mu;
        zetold =  zet;
        sold   =   s;
    %
        itto = 0;
        resinew = 2*residunorm;
        while resinew > residunorm && itto < 50
        itto = itto+1;
        x   =   xold + steg*dx;
        y   =   yold + steg*dy;
        z   =   zold + steg*dz;
        lam = lamold + steg*dlam;
        xsi = xsiold + steg*dxsi;
        eta = etaold + steg*deta;
        mu  = muold  + steg*dmu;
        zet = zetold + steg*dzet;
        s   =   sold + steg*ds;
        ux1 = upp-x;
        xl1 = x-low;
        ux2 = ux1.*ux1;
        xl2 = xl1.*xl1;
        uxinv1 = een./ux1;
        xlinv1 = een./xl1;
        plam = p0 + P'*lam ;
        qlam = q0 + Q'*lam ;
        gvec = P*uxinv1 + Q*xlinv1;
        dpsidx = plam./ux2 - qlam./xl2 ;
        rex = dpsidx - xsi + eta;
        rey = c + d.*y - mu - lam;
        rez = a0 - zet - a'*lam;
        relam = gvec - a*z - y + s - b;
        rexsi = xsi.*(x-alfa) - epsvecn;
        reeta = eta.*(beta-x) - epsvecn;
        remu = mu.*y - epsvecm;
        rezet = zet*z - epsi;
        res = lam.*s - epsvecm;
        residu1 = [rex' rey' rez]';
        residu2 = [relam' rexsi' reeta' remu' rezet res']';
        residu = [residu1' residu2']';
        resinew = sqrt(residu'*residu);
        steg = steg/2;
        end
      residunorm=resinew;
      residumax = max(abs(residu));
      steg = 2*steg;
      end
      if ittt > 198
        epsi;
        ittt;
      end
    epsi = 0.1*epsi;
    end
    xmma   =   x;
    ymma   =   y;
    zmma   =   z;
    lamma =  lam;
    xsimma =  xsi;
    etamma =  eta;
    mumma  =  mu;
    zetmma =  zet;
    smma   =   s;

end


%---------------------------------------------------------------------
%  This is the file kktcheck.m
%  Version Dec 2006.
%  Krister Svanberg <krille@math.kth.se>
%
function[residu,residunorm,residumax] = ...
kktcheck(m,n,x,y,z,lam,xsi,eta,mu,zet,s, ...
         xmin,xmax,df0dx,fval,dfdx,a0,a,c,d)
%
%  The left hand sides of the KKT conditions for the following
%  nonlinear programming problem are calculated.
%         
%      Minimize  f_0(x) + a_0*z + sum( c_i*y_i + 0.5*d_i*(y_i)^2 )
%    subject to  f_i(x) - a_i*z - y_i <= 0,  i = 1,...,m
%                xmax_j <= x_j <= xmin_j,    j = 1,...,n
%                z >= 0,   y_i >= 0,         i = 1,...,m
%*** INPUT:
%
%   m    = The number of general constraints.
%   n    = The number of variables x_j.
%   x    = Current values of the n variables x_j.
%   y    = Current values of the m variables y_i.
%   z    = Current value of the single variable z.
%  lam   = Lagrange multipliers for the m general constraints.
%  xsi   = Lagrange multipliers for the n constraints xmin_j - x_j <= 0.
%  eta   = Lagrange multipliers for the n constraints x_j - xmax_j <= 0.
%   mu   = Lagrange multipliers for the m constraints -y_i <= 0.
%  zet   = Lagrange multiplier for the single constraint -z <= 0.
%   s    = Slack variables for the m general constraints.
%  xmin  = Lower bounds for the variables x_j.
%  xmax  = Upper bounds for the variables x_j.
%  df0dx = Vector with the derivatives of the objective function f_0
%          with respect to the variables x_j, calculated at x.
%  fval  = Vector with the values of the constraint functions f_i,
%          calculated at x.
%  dfdx  = (m x n)-matrix with the derivatives of the constraint functions
%          f_i with respect to the variables x_j, calculated at x.
%          dfdx(i,j) = the derivative of f_i with respect to x_j.
%   a0   = The constants a_0 in the term a_0*z.
%   a    = Vector with the constants a_i in the terms a_i*z.
%   c    = Vector with the constants c_i in the terms c_i*y_i.
%   d    = Vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%     
%*** OUTPUT:
%
% residu     = the residual vector for the KKT conditions.
% residunorm = sqrt(residu'*residu).
% residumax  = max(abs(residu)).
%

    rex   = df0dx + dfdx'*lam - xsi + eta;
    rey   = c + d.*y - mu - lam;
    rez   = a0 - zet - a'*lam;
    relam = fval - a*z - y + s;
    rexsi = xsi.*(x-xmin);
    reeta = eta.*(xmax-x);
    remu  = mu.*y;
    rezet = zet*z;
    res   = lam.*s;
    %
    residu1 = [rex' rey' rez]';
    residu2 = [relam' rexsi' reeta' remu' rezet res']';
    residu = [residu1' residu2']';
    residunorm = sqrt(residu'*residu);
    residumax = max(abs(residu));

end


%----------------------------------------------------------
% PLOTTING FUNCTIONS
%----------------------------------------------------------

function plotwh(rho,ns,bc)
%PLOTWH     Plot the current state of the wheel and it spokes
% PATRICK KELLEY - UIUC - 2021

    [Xi,a] = physvar2Xi(rho,ns,bc);
    X = Xi(1:ns) + a(1);
    Y = Xi(ns+1:end) + a(2);
    plot([zeros(ns,1),X]',[zeros(ns,1),Y]','-k');
    hold on; axis equal; axis off;
    if ns<=40
        text(X,Y,num2str((1:ns)'));
    end
    plot([X;X(1)],[Y;Y(1)],'-k','LineWidth',3);
    % Show the axle location
    plot(0,0,'ro','MarkerSize',8,'MarkerFaceColor','r');
end


function plotopt(arc_c,A_T,L_T)
%PLOTOPT   Plots the optimal wheel shape given the
%           terrain parameters
% PATRICK KELLEY - UIUC - 2021

    if A_T == 0
        R_opt = arc_c/(pi);
        Q = 0:0.01:2*pi;
        X = R_opt.*cos(Q);
        Y = R_opt.*sin(Q);
    else
        d = 0.1*(L_T-1);
        dx = L_T/100;
        x = 0:dx:L_T;
        RR = d + A_T - (A_T/2).*(1-cos((2*pi/L_T).*x));
        TT = atan(sqrt(d).*tan(pi/L_T.*x)./sqrt(d+A_T))./(pi*sqrt(d)*sqrt(d+A_T)/L_T);
        t0 = -pi/2;
        X = RR.*cos(TT+t0);
        Y = RR.*sin(TT+t0);
    end
    plot(X,Y,'--g');
    hold on; axis equal; axis off;
end


function [xin] = continuation(xin,iter,inc,itersep,xmax)
%CONTINUATION  Function for handling a continuation method.
%               Takes xin as the parameter under continuation,
%               with iter being the current iteration in the 
%               loop. Inc is the increment which is added to
%               the parameter xin whenever iter is a multiple
%               of iterstep, provided xin < xmax.
%
% PATRICK KELLEY - UIUC - 2021

    if (mod(iter,itersep) == 0) && (xin < xmax)
        xin = xin + inc;
    end
end


function makegif(anim,fname,t)
%MAKEGIF  Turn an animation into a gif file.
%  INPUTS:
%   anim....the animation structure.
%   fname...the output gif's filename as a string, WITHOUT
%            a file extension (it is automagically appended).
%   t.......[Optional] framerate parameter. Default is 0, as
%             fast as possible.
%  OUTPUT:
%    The output gif is saved into the current file folder location.
%
% PATRICK KELLEY - UIUC - 2020

    if nargin < 3
        t = 0;
    end

    fname = [fname, '.gif'];
    framerate = t/length(anim);

    for j = 1:length(anim)
        im = frame2im(anim(j));
        [imind,cm] = rgb2ind(im,256);
        if j == 1
            imwrite(imind,cm,fname,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,fname,'gif','WriteMode','append','DelayTime',framerate);
        end
    end
end


function makepng(pic,fname)
%MAKEPNG  Turn a plot frame into a png file.
%  INPUTS:
%   pic.....the grabbed frame of the desired plot.
%   fname...the output png's filename as a string, WITHOUT
%            a file extension (it is automagically appended).
%  OUTPUT:
%    The output png is saved into the current file folder location.
%
% PATRICK KELLEY - UIUC - 2021

    fname = [fname, '.png'];
    im = frame2im(pic);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,fname,'png');
    
end
