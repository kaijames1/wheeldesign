clear; clc; close all;

parpool('local',12);

%==========================================================================
%=== TOPOLOGY OPTIMIZATION PARAMETERS =====================================
%==========================================================================
    % Material parameters
        E = 1e9; % Young's modulus
        v = 0.3; % Poisson's ratio

    % Optimization parameters
        p_SIMP = 2;         % Initial stiffness penalty exponent
        HbetaSO = 10;       % Heaviside filtering parameter for SO
        HbetaTO = 500;      % Heaviside filtering parameter for TO
        Vmax = 0.95;        % Max (initial) volume fraction
        VOm = 0.23;         % Max allowable volume fraction
        fr = 3.0;           % Filter radius, [elements]
        maxoutit = 2000;     % Maximum optimizer iterations

    % Set domain parameters for elasticity problem
        Ex = 160;            % number of elements in X direction
        Ey = Ex;            % number of elements in Y direction
        lx = 6.0;           % domain length in X direction
        ly = lx;            % domain length in Y direction
        n_layer = 2;        % number of layers in the FEA domain
        lne = Ex*Ey;        % number of elements per layer
        ne = lne*n_layer;   % total number of elements/design variables

    % Contact stiffness parameters
        sig = 0.5*lx/Ex;    % variance (width) of the Gaussian distribution
        Kcont0 = 1e4*E;     % nominal stiffness of contact interaction
        
    % Create default element
        scale = (lx/Ex + ly/Ey)/2; % Compute element length (assume square)
        Xe = scale*[0, 1, 0, 1];
        Ye = scale*[0, 0, 1, 1];
        k0 = Ke2Diso(E, v, Xe, Ye);
        
    % Initial material population
        rho = HeavisideInv(Vmax,HbetaTO)*ones(ne, 1);

    % Initialize boundary conditions
        [F_ext, constraints, ndof, axle] = define2D_multi(Ex, Ey);
        [nodeMap, Xm, Ym] = createMesh_multi(Ex, Ey, lx, ly);

    % Element connectivity filtering
        W_filt = density_filter(Ex, Ey, fr, Xm, Ym, nodeMap);
        rho1 = W_filt*rho(1:lne);       % layerwise filtered densities
        rho2 = W_filt*rho(lne+1:ne);
        rho_filt = [rho1;rho2];

    % Perform heaviside filtering
        [rho_fH,drHdr] = Heaviside(rho_filt,HbetaTO);

%% ========================================================================
% == SHAPE OPTIMIZATION PARAMETERS ========================================
%==========================================================================   
    %** ns : total number of spokes in the shape
    ns = 48;
        %-- Dependencies: DO NOT ALTER! ---------------------------
        %----nv : total number of design variables; the spokes plus
        %          all of the angles between them 
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
        wedgeang = 0.2*dQ;

      %** RADIAL CONTROL ******************************************
      %**** R_low : lower bound on the spoke length, [m]          *
      %**** R_upp : upper bound on the spoke length, [m]          *
      %************************************************************
        R_low = 0.3*lx;
        R_upp = 0.45*lx;     % MUST be <= lx/2

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
        L_T = (0.8*lx)*pi;

      %** SIMULATION PARAMETERS ***********************************
      %****** Fax : force on the axle in the vertical direction,  *
      %              [N]. Give as a positive scalar.              *
      %************************************************************
        Fax = 1e0;

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

    % Summation matrix for centroid
    Summer = zeros(2,2*ns);
    Summer(1,1:ns) = 1;    Summer(2,ns+1:end) = 1;

    %-- Generate the terrain structure and the arc length.
    [TX,dTXdX,arc_c] = genTX(A_T,L_T);

    %-- Compute the perimeter constraints.
    MinP = 2*pi*R_low;
    MaxP = pi*lx;   % To keep spokes inside the mesh!
    
    %--------------------------------------------------------------
    %-- RANDOM SHAPE GENERATION & FIRST COMPUTATIONS
    %--------------------------------------------------------------
    % Randomly initialize design variables
    x0 = rand(nv,1);

    tic;
    % Convert to physical variables
    xC = Wb*x0 + bb;

    % Compute the work done
    [W,dWdx] = rollsf(x0,ns,nv,Summer,Wb,bb,bc,Fax,TX,dTXdX,A_T,L_T,0);

    % Compute the initial perimeter
    [P,dPdrho] = chperim(xC,ns);
    % Perimeter sensitivities are known analytically if the phys.
    %  variables sent in were not convex-hulled
    dPdx = Wb'*dPdrho;

    % Generate pre-processable matrices
    kRmat = cell(ns,1);
    dRmatdx = cell(ns,1);
    Kcont = cell(ns,1);
    dKcdxRj = cell(ns,1);
    dKcdxQj = cell(ns,1);
    
    % Compute the absolute locations the spokes
    thetas = xC(ns+1:end) + bc;
    Wbvec = diag(Wb);
    
    parfor j = 1:ns
        %fprintf('Building Contact Matrix %d/%d...\n',j,ns);
        [kRmat{j},dRmatdtheta] = RblockGen(thetas(j),-pi/2);
        dRmatdx{j} = dRmatdtheta*Wbvec(ns+j);
        % Find the contact point and its sensitivities
        [XY,dXYdxR,dXYdxQ] = contactpt(j,ns,xC,Wb,bc);
        [Kcont{j},dKcdxRj{j},dKcdxQj{j}] = contactStiffnessNEW(XY, ...
                                    sig, Kcont0, ndof, lx, ly, Ex, Ey, ...
                                    dXYdxR,dXYdxQ);
        %clc;
    end
    
%%
% *************************************************************************
close; clc; tic;
fprintf('Iter.  KKT norm    f(x) = W   max(C)/Cmax   P/Pmin   V/Vmax  %c   %ct     %ct     %c\n',951,916,8747,946);
fprintf('----- ---------- ------------ ------------ -------- -------- - ------ ------ -----\n');



% INITIAL COMPLIANCE COMPUTATION
    C = zeros(ns, 1);

    dCdxtempS = zeros(nv,ns);
    dCdxtempT = zeros(ne,ns);
    
    drHdrB = drHdr(1:lne);
    drHdrW = drHdr(lne+1:ne);
    
    % Evaluate compliance and compliance sensitivities
    parfor i = 1:ns
        [C(i), dCdxT,dCdxS] = ...
            compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                                axle, Ex, Ey, nodeMap, k0, kRmat{i}, Kcont{i},...
                                ns, i, dRmatdx{i}, dKcdxRj{i}, dKcdxQj{i});
        dCdxtempT(:, i) = [((dCdxT(1:lne).*drHdrB)'*W_filt)';
                          ((dCdxT(lne+1:ne).*drHdrW)'*W_filt)'];              
        dCdxtempS(:,i) = dCdxS;
    end
    dCdx = [dCdxtempS;dCdxtempT];
    
    
    maxC = max(C);
    Cmax = 1.75*maxC; % Maximum allowable compliance
    
    % Volume constraint evaluation
    L = ones(lne, 1);
    Vol = sum(rho_fH);
    Vtot = ne;
    dVol(1:lne) = ((L.*drHdrB)'*W_filt)';
    dVol(lne+1:ne) = ((L.*drHdrW)'*W_filt)';
    
    
    
    
%% *************************************************************************
% Initialize parameters used in MMA problem
ntv = ne + nv;          % total no. of design vars
xval = [x0;rho];
m = ns + 3;        % ns compliance constraints + 2 perimeter cons + vol cons
epsimin = 1e-7;
xold1   = xval;
xold2   = xval;
xmin    = 0*ones(ntv, 1);
xmax    = ones(ntv, 1);
low     = xmin;
upp     = xmax;
c_mma   = 1000*ones(m, 1);
d_mma   = 1*ones(m, 1);
a0_mma  = 1;
a_mma   = 0*ones(m, 1);
outeriter = 0;
kkttol  = 1e-2;
kktnorm = kkttol + 1;


f0val = W;
df0dx = [dWdx;zeros(ne,1)];
fval = [MinP - P; P - MaxP; (C - Cmax)./Cmax; Vol/Vtot - VOm];
dfdx = [-[dPdx;zeros(ne,1)], [dPdx;zeros(ne,1)],...
          dCdx./Cmax, [zeros(nv,1); dVol'./Vtot]]; 
    

T = toc;
fprintf('%5d %10.4f %12.6f %12.6f %8.3f %8.3f %1d %6.2f %6.2f %5d\n',...
    outeriter,kktnorm,W,maxC/Cmax,P/MinP,Vol/Vtot,p_SIMP,T,T/60,HbetaSO);    

x_hist = zeros(length(xval),maxoutit);
C_hist = zeros(maxoutit,1);
V_hist = zeros(maxoutit,1);

tic;
%%
while outeriter < maxoutit && kktnorm > kkttol

    outeriter = outeriter+1;
    
    
    p_SIMP = continuation(p_SIMP,outeriter,1,50,5);
    HbetaSO = continuation(HbetaSO,outeriter,HbetaSO,30,500);

    
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(m,ntv,outeriter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval,dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma,HbetaSO,HbetaTO,nv);
    
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;

    % Grab the two different sets of design variables
    xvalW = xval(1:nv);         % Work design variables
    xvalC = xval(nv+1:end);     % Compliance design variables
    
    xC = Wb*xvalW + bb;
    % Compute the work done
    [W,dWdx] = rollsf(xvalW,ns,nv,Summer,Wb,bb,bc,Fax,TX,dTXdX,A_T,L_T,0);
    % Compute the initial perimeter
    [P,dPdrho] = chperim(xC,ns);
    dPdx = Wb'*dPdrho;
    
    thetas = xC(ns+1:end) + bc;
    
    parfor j = 1:ns
        [kRmat{j},dRmatdtheta] = RblockGen(thetas(j),-pi/2);
        dRmatdx{j} = dRmatdtheta*Wbvec(ns+j);
        % Find the contact point and its sensitivities
        [XY,dXYdxR,dXYdxQ] = contactpt(j,ns,xC,Wb,bc);
        [Kcont{j},dKcdxRj{j},dKcdxQj{j}] = contactStiffnessNEW(XY, ...
                                    sig, Kcont0, ndof, lx, ly, Ex, Ey, ...
                                    dXYdxR,dXYdxQ);
    end
    
    % Filter the densities layer-wise
    rho1 = W_filt*xvalC(1:lne);
    rho2 = W_filt*xvalC(lne+1:ne);
    % Perform the Heaviside filtering
    [rho_fH,drHdr] = Heaviside([rho1;rho2],HbetaTO);
    % Grab the Heaviside sensitivities for each respective layer
    drHdrB = drHdr(1:lne);
    drHdrW = drHdr(lne+1:ne);

    % Evaluate compliance and compliance sensitivities
    parfor i = 1:ns
        [C(i), dCdxT,dCdxS] = ...
            compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                                axle, Ex, Ey, nodeMap, k0, kRmat{i}, Kcont{i},...
                                ns, i, dRmatdx{i}, dKcdxRj{i}, dKcdxQj{i});
        dCdxtempT(:, i) = [((dCdxT(1:lne).*drHdrB)'*W_filt)';
                          ((dCdxT(lne+1:ne).*drHdrW)'*W_filt)'];              
        dCdxtempS(:,i) = dCdxS;
    end
    dCdx = [dCdxtempS;dCdxtempT];
    
    % Evaluate volume and sensitivities
    Vol = sum(rho_fH);   
    dVol(1:lne) = ((L.*drHdrB)'*W_filt)';
    dVol(lne+1:ne) = ((L.*drHdrW)'*W_filt)';
    
    % Objective and constraint functions, sensitivities updated
    f0val = W;
    df0dx = [dWdx;zeros(ne,1)];
    fval = [MinP - P; P - MaxP; (C - Cmax)./Cmax; Vol/Vtot - VOm];
    dfdx = [-[dPdx;zeros(ne,1)], [dPdx;zeros(ne,1)],...
              dCdx./Cmax, [zeros(nv,1); dVol']./Vtot];    
    
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = kktcheck(m,ntv,xmma,ymma,zmma,...
        lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx',...
        a0_mma,a_mma,c_mma,d_mma);
    
    %%%% Output information for current iteration
    T = T + toc;
    fprintf('%5d %10.4f %12.6f %12.6f %8.3f %8.3f %1d %6.2f %6.2f %5d\n',...
        outeriter,kktnorm,W,max(C)/Cmax,P/MinP,Vol/Vtot,p_SIMP,toc,T/60,HbetaSO); 
    tic;

    % Plot component topologies
    if max(outeriter == [100, 200, 400, 800, 1200, 1600, maxoutit])
        outname = ['MDO_iter',num2str(outeriter),'_'];

        figure(2);
            structplot(rho_fH(1:lne),Ex,Ey,lx,ly);
            title(['Barrow Topology | Iteration: ',num2str(outeriter)]);
            makepng(getframe(2),[outname,'Barrow']);

        figure(3);
            structplot(rho_fH(lne+1:end),Ex,Ey,lx,ly);
            title(['Wheel Topology | Iteration: ',num2str(outeriter)]);
            makepng(getframe(3),[outname,'Wheel']);

        figure(4);
            plotwh(xC,ns,bc,lx,ly); hold on;
            plotopt(MinP,A_T,L_T,lx); plotopt(MaxP,A_T,L_T,lx); hold off;
            title(['Wheel Shape | Iteration: ',num2str(outeriter)]);
            makepng(getframe(4),[outname,'Shape']);
    end
    
    figure(5)
        subplot(1,3,1)
            structplot(rho_fH(1:lne),Ex,Ey,lx,ly);
        subplot(1,3,2)
            structplot(rho_fH(lne+1:end),Ex,Ey,lx,ly);
        subplot(1,3,3)
            plotwh(xC,ns,bc,lx,ly); hold on;
            plotopt(MinP,A_T,L_T,lx); plotopt(MaxP,A_T,L_T,lx); hold off;
        title(['MDO | Iteration: ',num2str(outeriter)]);
    allanim(outeriter) = getframe(5);
    
    % Store variable histories
    x_hist(:, outeriter) = xval;
    C_hist(outeriter) = max(C);
    V_hist(outeriter) = Vol;

end

%%

makegif(allanim,'MDO_GIF');
makepng(allanim(end),'MDO_FINAL');

structplot(rho_fH(1:lne),Ex,Ey,lx,ly,'MDO_FinalBarrow.png');
structplot(rho_fH(lne+1:end),Ex,Ey,lx,ly,'MDO_FinalWheel.png');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -----  |   |  \   |  /----  _____  -----  /---\  \   |  /----   %
%   |      |   |  |\  |  |        |      |    |   |  |\  |  |       %
%   ----   |   |  | \ |  |        |      |    |   |  | \ |  \---\   %
%   |      |   |  |  \|  |        |      |    |   |  |  \|      |   %
%   |      \___/  |   \  \____    |    -----  \___/  |   \  ____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Rblock,dRblockdQ] = RblockGen(Q,t0)
%RBLOCKGEN   Generates the 8x8 matrix for rotating an element's
%             stiffness matrix by angle Q+t0 WRT to the global
%             reference frame. Can return sensitivities of the
%             matrix WRT Q if so desired, as the second output.
%
% PATRICK KELLEY - UIUC - 2021

    if nargin < 2; t0 = 0; end

    Rmat = [cos(Q+t0), -sin(Q+t0); sin(Q+t0), cos(Q+t0)];
    Rblock = zeros(8, 8);
    
    if nargout > 1
        dRmat = [-sin(Q+t0),  -cos(Q+t0); cos(Q+t0), -sin(Q+t0)];
        dRblockdQ = zeros(8, 8);
    end
    
    for i = 1:4
        Rblock(2*i-1:2*i, 2*i-1:2*i) = Rmat;
        if nargout > 1
            dRblockdQ(2*i-1:2*i, 2*i-1:2*i) = dRmat;
        end
    end    
end


function [Z,dZdR,dZdQ] = contactpt(j,ns,xC,Wb,bc)

    % Contact point X and Y coordinate
    Z    = [xC(j)*cos(xC(j+ns) + bc(j));
            xC(j)*sin(xC(j+ns) + bc(j))];
    % Sensitivity of contact pt WRT radial design variable
    dZdR = [ Wb(j,j)*cos(xC(j+ns) + bc(j));
             Wb(j,j)*sin(xC(j+ns) + bc(j))];
     % Sensitivity of contact pt WRT angular design variable
    dZdQ = [-xC(j)*sin(xC(j+ns) + bc(j))*Wb(j+ns,j+ns);
             xC(j)*cos(xC(j+ns) + bc(j))*Wb(j+ns,j+ns)];

end


function [Kcont,dKcR,dKcQ] = contactStiffnessNEW(pin, sig, K0, ndof, ...
                                lx, ly, Ex, Ey, dpindxR, dpindxQ)
% Function computes the stiffness contributions due to the wheel's contact
% with the ground.
%
%  INPUTS:
%   pin.......[X;Y] location of the current contact point, origin at (0,0)
%   sig.......Spatial width of the gaussian distribution
%   K0........Scalar stiffness multiplier
%   ndof......Number of Degrees of Freedom in the mesh
%   lx, ly....Length of the domain in the X, Y directions
%   Ex, Ey....Number of elements in the X, Y directions
%   dpindxR...Sens of pin location WRT xR, (2 x 1)
%   dpindxQ...Sens of pin location WRT xQ, (2 x 1)
%
%  OUTPUTS:
%   Kcont.....Contact stiffness matrix for wheel layer
%   dKcR......Sensitivities matrix WRT xR design variable
%   dKcQ......Sensitivities matrix WRT xQ design variable
%
% PATRICK KELLEY - UIUC - 2021

    % Initialize modified global (layer) stiffness and sens. matrices
    Kcont = zeros(ndof,1);  
    dKcR = zeros(ndof,1);
    dKcQ = zeros(ndof,1);

    
    % Want to see the distribution of the contact stiffness? Toggle to 1
    shoplot = 0;
    
    % Spatial width of the elements in each direction
    dx = lx/Ex; dy = ly/Ey;

    kl = K0*ones(2,1);             % Link stiffness
    coeff = sqrt(2*pi)*sig;     % Coefficient used in Gaussian function evaluation
    
    % Array used for visualizing placement of stiffness, if desired
    if shoplot == 1; vis = zeros(Ey+1, Ex+1); end
    
    for j = 1:Ey+1
        
        Y = -ly/2 + (j-1)*dy;       % Y-location of current node
        
        for k = 1:Ex+1
         
            X = -lx/2 + (k-1)*dx;   % X-location of current node
            
            % Squared distance between FEM node and spoke tip
            dist2 = (X-pin(1))^2 + (Y-pin(2))^2;
            
            % Gaussian distribution for stiffness value
            Phi = coeff*exp(-dist2/(2*sig^2));
            k_i = Phi*kl;

            % Store for plotting in contact map, if desired
            if shoplot == 1; vis(j, k) = vis(j, k) + Phi; end
            
            % Find where the stiffness is to be added to the matrix
            nodeNum = (j-1)*(Ex+1) + k;
            edofs = [2*nodeNum-1, 2*nodeNum];
            % Add it in!
            Kcont(edofs,1) = Kcont(edofs,1) + k_i;
            
            % SENSITIVITIES of the contact matrix WRT shape design vars
                dKcR(edofs,1) = dKcR(edofs,1) + ...
                            k_i*((X-pin(1))*dpindxR(1)+ ...
                                 (Y-pin(2))*dpindxR(2))/sig^2;
                dKcQ(edofs,1) = dKcQ(edofs,1) + ...
                            k_i*((X-pin(1))*dpindxQ(1)+ ...
                                 (Y-pin(2))*dpindxQ(2))/sig^2;
        end
    end
    
    % Sparsify the matrices for a bit better performance
    Kcont = spdiags(Kcont,0,ndof,ndof);
    dKcR = spdiags(dKcR,0,ndof,ndof);
    dKcQ = spdiags(dKcQ,0,ndof,ndof);

    % Plot contact map, if desired
    if shoplot == 1
        figure(6); cla; axis equal;
        [Xgrid, Ygrid] = meshgrid(-lx/2:dx:lx/2, -ly/2:dy:ly/2);
        surf(Xgrid, Ygrid, vis); hold on
        plot(pin(1), pin(2), 'm*'); hold on
        view([0 90]);
    end
end


function [C, dCT, dCS] = ...
    compliance2D_multi(rho, p_SIMP, F, cons, axle, Ex, Ey, nodeMap, k0, Rmat, Kcont,...
                        ns,Citer,dRmat,dKcdxR,dKcdxQ)
% function calculates the displacements and compliance of the topology 
% defined by the relative density matrix rho, subject to the applied force F 
% the vector 'constraints' contains the numbers of the constrained 
% degrees of freedom and
% dimx and dimy are the number of elements along the x and y dimensions 
% of the domain
%
% INPUTS:
%   rho.......(2nepl x 1) filtered density variable vector
%   p_SIMP....(scalar) SIMP penalty exponent
%   F.........(ndofpl x 1) forces on the barrow layer vector
%   cons......constraint location vector
%   axle......axle node location vector
%   Ex........Number of elements in the X
%   Ey........"                       " Y
%   nodeMap...Matrix for mapping nodes to elements
%   k0........Default stiffness matrix, (8x8)
%   Rmat......Current rotation matrix R(theta_j), (8x8)
%   Kcont.....Contact stiffness matrix, assembled, (nnpl x nnpl)
%   ns........Number of spokes in the shape
%   Citer.....Current number of the compliance computation, for computing
%              where sensitivities go for the contact coupling of design
%   dRmat.....Sensitivities of rotation matrix WRT spoke design var Q, (8x8)
%   dKcdxR....(nnpl x nnpl) Contact stiffness matrix sensitivities WRT the
%              design variable controlling the Citer-th spoke length
%   dKcdxQ....(nnpl x nnpl) Contact stiffness matrix sensitivities WRT the
%              design variable controlling the Citer-th angular location
%
%  OUTPUTS:
%   C.........Compliance value with Citer-th spoke in contact with ground
%   dC........Sensitivities of said compliance WRT design variables
%
% KAI JAMES | PATRICK KELLEY - UIUC - 2021


    ndof = (Ey+1)*(Ex+1)*2;
    el_dim = 2;

    lne = Ex*Ey; % number of elements per layer
    ne = 2*lne;  % total number of elements

    % Rotated element stiffness matrix
    kR = Rmat*k0*Rmat';
    [K1, K2] = global_stiffness2D_multi(rho, p_SIMP, kR, k0, el_dim, Ex, Ey, nodeMap);  

    % Define constraints/supports
    ndofGlob = 2*ndof-2; % # of dofs in global system with axle constraint
    freedofs = setdiff(1:ndofGlob, cons);

    Uglob = zeros(ndofGlob, 1);
    Fglob = sparse(ndofGlob, 1);

    % Assemble global stiffness matrix from layer stiffness matrices
    L1dofs = 1:ndof; % layer 1 dofs
    L2dofs = [ndof+1:ndof+axle(1)-1, axle, ndof+axle(1):ndofGlob]; % all layer 2 dofs expressed in global indexing
    Kglob = sparse(ndofGlob, ndofGlob);
    Kglob(L1dofs, L1dofs) = Kglob(L1dofs, L1dofs) + K1;
    Kglob(L2dofs, L2dofs) = Kglob(L2dofs, L2dofs) + K2 + Kcont; 

    % Solve linear systems and re-assemble global matrices/vectors
    Fglob(L1dofs) = F;
    Uglob(freedofs) = Kglob(freedofs, freedofs)\Fglob(freedofs, :);
    C = Fglob'*Uglob;

    % Extract layer-specific displacement vectors
    U1 = Uglob(L1dofs);
    U2 = Uglob(L2dofs);

    % Sens of C WRT shape design vars x and then topology design vars rho
    dCT = zeros(ne,1);
    dCS = zeros(ns*2,1);
    edof = zeros(2*el_dim^2, 1);
    
    % COMPUTE COMPLIANCE SENSITIVITIES
    for i = 1:lne

        % Global degrees of freedom affected by element (i, j) 
        for m = 1:el_dim^2
            edof(2*m - 1) = 2*nodeMap(i,m) - 1;
            edof(2*m) = 2*nodeMap(i,m);
        end

        % SENSITIVITY WRT filtered element denisity design var, first for
        %  the barrow layer, then the wheel layer
        dCT(i) = -p_SIMP*rho(i)^(p_SIMP-1)*U1(edof)'*k0*U1(edof);
        dCT(i+lne) = -p_SIMP*rho(i+lne)^(p_SIMP-1)*U2(edof)'*kR*U2(edof);
        
        % SENSITIVITY WRT spoke length and contact point
        dCS(Citer) = dCS(Citer) - U2(edof)'*dKcdxR(edof,edof)*U2(edof);
        
        % SENSITIVITY WRT angular spoke location and contact point
        dCS(Citer+ns) = dCS(Citer+ns) - U2(edof)'*(...
            dRmat*(k0)*Rmat' + Rmat*(k0)*dRmat' + ...
            dKcdxQ(edof,edof))*U2(edof)*rho(i)^(p_SIMP);
    end
end


function [nodeMap, X, Y] = createMesh_multi(Ex, Ey, lx, ly)
%CREATEMESH_MULTI  Create the node mapping of a single Ex x Ey
%                   layer (assuming both layers have identical
%                   numbering). Also assumes linear elements.
%
%
% KAI JAMES - PATRICK KELLEY | UIUC | 2020

    nel = Ex*Ey;                % Number of elements per layer
    nodeMap = zeros(2*nel, 4);  % Preallocate

    % Build the node connectivity matrix
    for j = 1:Ey
        for i = 1:Ex
            ind = (j-1)*Ex + i;
            nodeMap(ind, 3) = (j-1)*(Ex+1) + i;
            nodeMap(ind, 4) = (j-1)*(Ex+1) + i+1;
            nodeMap(ind, 1) = (j)*(Ex+1) + i;
            nodeMap(ind, 2) = (j)*(Ex+1) + i+1;
        end
    end

    % Create vectors with the X,Y locations of all the nodes
    if nargout > 1
        dx = lx/Ex;
        dy = ly/Ey;
        X = zeros((Ex+1)*(Ey+1),1);
        Y = zeros((Ex+1)*(Ey+1),1);
        for j = 1:Ey+1
            for i = 1:Ex+1
                ind = (j-1)*(Ex+1) + i;
                X(ind, 1) = (i-1)*dx;
                Y(ind, 1) = ly - (j-1)*dy;
            end
        end
    end
end


function [F_ext, constraints, ndof, axle] = define2D_multi(Ex, Ey)

% Function conputes the consistent force vector used for external forces
% acting on layer 1 ONLY
% Function also identifies the constrained dofs in layer 1

    ndof = 2*(Ex+1)*(Ey+1); % number of dof's (per layer)
    F_ext = zeros(ndof, 1);


    % Force due to evenly distributed dead load (layer 1)

    tw = 1e5; % total weight of dead load

    for i = 1:Ex/2
        F_ext(Ey*(Ex+1)*2 + 2*i) = F_ext(Ey*(Ex+1)*2 + 2*i) - tw/Ex/2;
        F_ext(Ey*(Ex+1)*2 + 2*(i+1)) = F_ext(Ey*(Ex+1)*2 + 2*(i+1)) - tw/Ex/2;
    end

    F_ext = sparse(F_ext);

    % Vertically constraint degrees of freedom at handle location (assume user
    % will hold the handle vertically steady and this vertical component of
    % the pulling force will perform negligible work

    % Vertically and horizontally constrain the dof where handle is connected
    constraints(1) = 2*(Ex+1)*(Ey+1) - 1;
    constraints(2) = 2*(Ex+1)*(Ey+1);

    % Identify dof's corresponding to axle constraint
    % Axle is currently placed at the center of the layer, but could be moved
    % back
    ann = (Ey/2)*(Ex+1) + Ex/2 +1; % "axle node number"
    axle = [2*ann-1, 2*ann];

end


function W = density_filter(Ex, Ey, nh, X, Y, nodeMap)
% function generates a coefficient matrix used for filtering of element
% densities
% W contains the weight coefficients for a linear filter with radius given
% by the parameter 'rad'
% compute element centroid locations

    n_e = Ex*Ey;
    X_cent = zeros(n_e, 1);
    Y_cent = zeros(n_e, 1);

    for i = 1:n_e

        X_cent(i) = sum(X(nodeMap(i, :)))/4;
        Y_cent(i) = sum(Y(nodeMap(i, :)))/4;

    end    

    % find 'neighborhood' size (i.e. radius - measured in elements - of the filter neighborhood)
    % based on characteristic lenght of elements
    scale = max(X_cent(2)-X_cent(1), Y_cent(2)-Y_cent(1));
    rad = nh*scale;

    W = zeros(n_e, n_e);

    for i = 1:n_e

        weight_sum = 0;

        for iw = 1:n_e

            dist = sqrt((X_cent(i) - X_cent(iw))^2 + (Y_cent(i) - Y_cent(iw))^2);
            coeff = rad - dist;

            if coeff > 0
                W(i, iw) = coeff;
                weight_sum = weight_sum + coeff;
            end

        end

        % normalize weights
        W(i, :) = W(i, :)/weight_sum;

    end

    W = sparse(W);

end


function [K1, K2] = global_stiffness2D_multi(rho, p, kR, k0, el_dim, Ex, Ey, nodeMap)
% Script used to generate a global stiffness matrix
% function assembles the global stiffness matrices for each layers for the 
% topology defined by the relative density matrix `layout' whose elements 
% have specific stiffness matrices as given in the 4 dimensional array 
% k_elem4
% Note: el_dim = dimension of elements measured in nodes across
%
%*** NOTE: ASSUME 2D QUAD ELEMENTS ***

    ne = Ex*Ey;         % Number of elements in the mesh
    rho_min = 1e-6;     % Minimum stiffness value, to prevent singularness

    % Reshape the stiffness matrix into a vector for more efficient
    %  storage in creating the stiffness matrix
        k1v = reshape(k0,64,1);
        k2v = reshape(kR,64,1);

    % PREALLOCATIONS
    I = zeros(64,ne);       % Row indices storage vector
    J = zeros(64,ne);       % Column indices storage vector
    K1vec = zeros(64,ne);   % Layer 1 (Carriage) stiffness storage matrix
    K2vec = zeros(64,ne);   % Layer 2 (Wheel) stiffness storage matrix
    edof = zeros(8, 1);     % Mapping of element DOF vector
    
    for i = 1:ne
        % Compute where the element stiffness should be stored in the
        %  layer's global positioning
        for m = 1:el_dim^2
            edof(2*m - 1) = 2*nodeMap(i,m) - 1;
            edof(2*m) = 2*nodeMap(i,m);
        end
        % Store the appropriate indexing
        I(:,i) = reshape(repmat(edof,1,8)',64,1);
        J(:,i) = repmat(edof,8,1);
        % Store the element stiffnesses in the abbreviated matrix
        K1vec(:,i) = K1vec(:,i) + k1v*((rho(i))^p + rho_min);
        K2vec(:,i) = K2vec(:,i) + k2v*((rho(ne+i))^p + rho_min);
    end
    
    
    % Create the sparse matrix using the indexing and non-zero stiffness
    %  element vectors
    K1 = sparse(I,J,K1vec);
    K2 = sparse(I,J,K2vec);
end


function [U,dUdx] = Heaviside(x,Hbeta)
%HEAVISIDE  Smooth Heaviside function as used for
%            topology optimization.
%
%  INPUTS:
%   x.......Value of radius-filtered element denisities.
%   Hbeta...Smoothing exponent. HEAVISIDE approaches
%            discontinuous step function as Hbeta -> inf.
%  OUTPUT:
%   U.......Filtered element densities.
%   dUdx....[Optional] Sensitivity of HEAVISIDE function
%            to the input densities.
%
% PATRICK KELLEY - UIUC - 2020

    U = 1 - exp(-x*Hbeta) + x*exp(-Hbeta);
    
    if nargout > 1
        dUdx = exp(-x*Hbeta)*Hbeta + exp(-Hbeta);
    end
end


function Hi = HeavisideInv(v0,beta)


    f = @(x) exp((1-x)*beta) - x - exp(beta)*(1-v0);
    Hi = fzero(f,1/beta);
end


function [ke, vol] = Ke2Diso(E, v, X, Y)
% Function generates the stiffness matrix for an arbitrary quadrilateral 
% element given the Young's modulus E, Poisson's ratio v, and vectors
% containing the X and Y coordinates of the nodes of the element
% The function uses an isoparametric formulation

    ke = zeros(8, 8);
    vol = 0;
    D = constitutive(E, v);

    for j = 1:2
        for i = 1:2
            imap = (2*i-3)/sqrt(3);
            jmap = (2*j-3)/sqrt(3);
            [B, detJ] = StrainDisp(X, Y, imap, jmap);
            ke = ke + B'*D*B*detJ;
        end
    end
end


function dN = shapeFunctions(s, t)
% Shape Functions
% This function evaluates the DERIVATIVES of the shape function
% given the natural coordinates, s, t

    for j = 1:2
        for i = 1:2
            ind = (j-1)*2 + i;
            imap = 2*i-3;
            jmap = 2*j-3;
            dN(ind, 1) = imap*(1+jmap*t)/4;
            dN(ind, 2) = (1+imap*s)*jmap/4;
            
        end
    end
end


function [B, detJ] = StrainDisp(X_el, Y_el, s, t)
% Strain-Displacement Matrix
    dN = shapeFunctions(s, t);

    dx = [0, 0];
    dy = [0, 0];

    for i = 1:4
        dx = dx + dN(i, :)*X_el(i);
        dy = dy + dN(i, :)*Y_el(i);
    end

    J = [dx', dy'];
    detJ = det(J);

    B = zeros(3,8);
    
    for i = 1:4
        sp = (i-1)*2;

        M1 = [dN(i, 1), dy(1); dN(i, 2), dy(2)];
        M2 = [dx(1), dN(i, 1); dx(2), dN(i, 2)];

        B(1, sp+1) = det(M1);
        B(2, sp+2) = det(M2);

        B(3, sp+1) = det(M2);
        B(3, sp+2) = det(M1);
    end

    B = B/detJ;

end


function [D, M] = constitutive(E, v)
% Constitutive Matrix 
  % Generate constitutive matrix D given Young's Modulus E, and
  % Poisson's ratio v
  
  D(1, 1) = E/(1-v^2);
  D(1, 2) = v*E/(1-v^2);
  D(2, 1) = D(1, 2);
  D(2, 2) = D(1, 1);
  D(3, 3) = E/(2*(1+v));
  
  if nargout > 1
      M(1, 1) = 1;
      M(1, 2) = -1/2;
      M(2, 1) = M(1, 2);
      M(2, 2) = M(1, 1);
      M(3, 3) = 3;
  end
  
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SO FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TX,dTXdX,arc_c] = genTX(A_T,L_T)
    if A_T == 0
        TX = @(X) 0;
        dTXdX = @(X) 0;
        arc_c = L_T/2;
    else
        TX = @(X) A_T.*sin(pi.*X./L_T).^2;
        dTXdX = @(X) 2*A_T*sin(pi*X/L_T)*cos(pi*X/L_T)*pi/L_T;
        dx = 1e-3;
        Gx = 0:dx:L_T/2;
        arc_c = trapz(sqrt(1+(A_T*pi/L_T*sin(2*pi/L_T*Gx)).^2))*dx;
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
%              present; the sweep the design variable x can 
%              control via (j-1)dQ - wedge/2 + wedge*x_j if
%              wedge <= dQ.
%  OUTPUTS:
%   W.........2ns x 2ns matrix
%   b.........2ns x 1 constant vector
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
        
    end
    
    
    % CREATE THE GLOBAL ARRAYS
    W = blkdiag(Wr,Wq);
    b = [br;bq];
end


function [Xi,a,dXidrho,dadx] = physvar2Xi(rho,ns,bc)
%PHYSVAR2Xi   Convert radial/angular physical variables
%             into complex numbers using polar form.
%
% INPUTS:
%   rho...(2ns x 1) physical variables. First ns are the
%           spoke lengths, next ns are converted angles [rad]
%   ns....number of spokes in the shape
%   bc....(ns x 1) vector for shifting the angles
%
% OUTPUTS:
%   Z.....(ns x 1) locations of the spoke tips in 
%           the complex plane
%   dZdrh.(ns x 2ns) sensitivity matrix
%
% PATRICK KELLEY - UIUC 2020

    Xi = zeros(2*ns,1);
    
    R = rho(1:ns);
    Q = rho(ns+1:end);
    
    a = [R(1)*cos(Q(1) + bc(1));R(1)*sin(Q(1) + bc(1))];
    
    Xi(1:ns)     = R.*cos(Q + bc) - R(1)*cos(Q(1) + bc(1));
    Xi(ns+1:end) = R.*sin(Q + bc) - R(1)*sin(Q(1) + bc(1));

    
    if nargout >= 3
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
    
    if nargout == 4
        dadx = zeros(2,2*ns);
        dadx(1,1) = cos(Q(1) + bc(1));
        dadx(1,ns+1) = -R(1)*sin(Q(1) + bc(1));
        dadx(2,1) = sin(Q(1) + bc(1));
        dadx(2,ns+1) = R(1)*cos(Q(1) + bc(1));
    end    
end


function plotwh(xC,ns,bc,lx,ly)
    % Convert to cartesian
    Ox = lx/2; Oy = ly/2;
    X = cos(xC(ns+1:end) + bc).*xC(1:ns) + Ox;
    Y = sin(xC(ns+1:end) + bc).*xC(1:ns) + Oy;
    % Plot
    plot([Ox+zeros(ns,1),X]',[Oy+zeros(ns,1),Y]','-ok');
    hold on; axis equal; axis off;
    plot([X;X(1)],[Y;Y(1)],'-ok');
    % Show the axle location
    plot(Ox,Oy,'ro','MarkerSize',8,'MarkerFaceColor','r');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MMA OPTIMIZER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%---------------------------------------------------------------------
end


%-------------------------------------------------------
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,HbetaSO,HbetaTO,nsv)
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
%--------------------------------------------------------------------------

    epsimin = 10^(-7);
    raa0 = 0.00001;
    albefa = 0.1;
    % Changed asyinit so that the two Heaviside for the
    %   two optimization regimes aren't dependent
    %     asyinit = 0.5/(1+H_beta);
        asyinit = ones(n,1);
        asyinit(1:nsv) = 0.5/(1+HbetaSO);
        asyinit(nsv+1:end) = 0.5/(1+HbetaTO);
        %
    asyincr = 1.02;      % Modified by PK, was 1.2
    asydecr = 0.98;      % "            ", was 0.7
    eeen = ones(n,1);
    eeem = ones(m,1);

    % Calculation of the asymptotes low and upp :
    if iter < 2.5
        low = xval - asyinit;
        upp = xval + asyinit;
    else
        zzz = (xval-xold1).*(xold1-xold2);
        factor = eeen;
        factor((zzz > 0)) = asyincr;
        factor((zzz < 0)) = asydecr;
        low = xval - factor.*(xold1 - low);
        upp = xval + factor.*(upp - xold1);
        lowmin = xval - 10*(xmax-xmin);
%         lowmax = xval - min(0.01,1/H_beta)*(xmax-xmin);     % Modified by PK
%         uppmin = xval + min(0.01,1/H_beta)*(xmax-xmin);     % "
            % Modified to decouple the heavisiding of the regimes
            lowmax = xval - 2*asyinit;     % Modified by PK
            uppmin = xval + 2*asyinit;     % "
            %
        uppmax = xval + 10*(xmax-xmin);
        low = max(low,lowmin);
        low = min(low,lowmax);
        upp = min(upp,uppmax);
        upp = max(upp,uppmin);
    end

    % Calculation of the bounds alfa and beta :
    zzz = low + albefa*(xval-low);
    alfa = max(zzz,xmin);
    zzz = upp - albefa*(upp-xval);
    beta = min(zzz,xmax);

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
    p0 = max(df0dx,0);
    q0 = max(-df0dx,0);
    pq0 = 0.001*(p0 + q0) + raa0*xmamiinv;
    p0 = p0 + pq0;
    q0 = q0 + pq0;
    p0 = p0.*ux2;
    q0 = q0.*xl2;
    %
    P = max(dfdx,0);
    Q = max(-dfdx,0);
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
%-------------------------------------------------------------
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PLOTTING AND MISCELLANEOUS FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function structplot(rho,nex,ney,Lx,Ly,fname)

    rrho = zeros(ney, nex);

    for j = 1:ney
        for i = 1:nex
            ind = (j-1)*nex + i;
            rrho(ney+1-j, i) = rho(ind);
        end
    end

    imagesc([0 Lx],[0 Ly],1-rrho);    
    colormap('gray');
    caxis([0 1]);
    grid off;
    view(0,90);
    xlim([-0.05 1.05]*Lx);
    ylim([-0.05 1.05]*Ly);
    axis equal;
    axis off;
    
    if nargin > 5
        imwrite(1-rrho,fname)
    end
end

function plotopt(arc_c,A_T,L_T,lx)
%PLOTOPT   Plots the optimal wheel shape given the
%           terrain parameters
% PATRICK KELLEY - UIUC - 2021

    if A_T == 0
        R_opt = arc_c/(2*pi);
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
    plot(X+lx/2,Y+lx/2,'--g');
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