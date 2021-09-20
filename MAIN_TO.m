% TOPOLOGY OPTIMIZATION MAIN PROGRAM
%  for rigid-body, left-to-right rolling "wheels" over various
%   terrain types. optimally synthesizing a wheelbarrow mechanism.  
%   Optimizes the wheel's internal structural topology, as well
%   as the topology of the carriage.
%  Objective function type: MASS MINIMIZATION,
%   i.e. f = sum(rho_fH)
%  subject to compliance constraints, F'd <= Cmax
%
% PATRICK KELLEY | KAI JAMES - UIUC - 2020
%
clear; clc; close all;

usepar = 0; % toggle whether or not to use parallel functions

if usepar; parpool('local',12); end
%
%**************************************************************************
%** USER INPUTS ***********************************************************
%**************************************************************************
%
%   TOPOLOGY OPTIMIZATION PARAMETERS
%==========================================================================
%**Material parameters, as used in a plane-stress FEM formulation
    E = 1e9;            % Young's modulus
    v = 0.3;            % Poisson's Ratio

%**Optimization parameters
    p_SIMP = 2;         % Initial stiffness penalty exponent
    maxSIMP = 5;        % Maximum value of SIMP exponent in continuation
    H_beta = 500;       % Heaviside filtering parameter (initially)
    Vmax = 0.99;        % Max (initial) volume fraction
    fr = 3.0;           % Filter radius, [elements]
    maxoutit = 1000;    % Maximum optimizer iterations
    kkttol  = 1.0e-3;   % Convergence tolerance

%**Domain parameters for elasticity problem
    Ex = 80;           % number of elements in X direction
    lx = 6.0;           % domain length in X direction
    n_layer = 2;        % number of layers in the FEA domain
    % DEPENDENCIES
        Ey = Ex;        % number of elements in Y direction
        ly = lx;        % domain length in Y direction
        lne = Ex*Ey;    % number of elements per layer
        ne = lne*n_layer; % total number of finite elements/design vars
    
%**Contact stiffness parameters
    sig = 0.5*lx/Ex;    % variance (width) of the Gaussian distribution
    Kcont0 = 1e4*E;     % nominal stiffness of contact interaction

%   SHAPE PARAMETERS
%==========================================================================
%**General wheel properties
    ns = 16; % number of sectors (spokes) in the wheel
    n = 2*ns; % number of design parameters
    dQ = 2*pi/ns; % Angle of each wheel sector

%**Set design variables for CIRCLE  (comment/uncomment as desired)
     xC = 2.7*ones(2*ns,1);
     xC(ns+1:end) = dQ;
%**Set design variables for ELLIPSE (uncomment/comment as desired) (ns=48)
%      [xC] = ellipsegen(ns);
  
     
     
%--------------------------------------------------------------------------
%   COMPUTE INITIAL DEPENDENCIES
%--------------------------------------------------------------------------
%**Create default element
    scale = (lx/Ex + ly/Ey)/2; % Compute element length (assume square)
    Xe = scale*[0, 1, 0, 1];
    Ye = scale*[0, 0, 1, 1];
    k0 = Ke2Diso(E, v, Xe, Ye);
%**Initialize boundary conditions
    [F_ext, constraints, ndof, axle] = define2D_multi(Ex, Ey);
    [nodeMap, Xm, Ym] = createMesh_multi(Ex, Ey, lx, ly);

%**Generate the linear density filtering matrix
    W_filt = density_filter(Ex, Ey, fr, Xm, Ym, nodeMap); 
    % W_filt = eye(lne); % #nofilter

%--------------------------------------------------------------------------
%   GENERATE PRE-PROCESSABLE MATRICES
%--------------------------------------------------------------------------        
%**The two matrix types that can be created prior to the optimization are
%   the 8x8 "R_block" matrix that rotates each element default stiffness
%   matrix in accordance to the orientation of the contacting spoke. Also,
%   the contacting point matrices can be found beforehand, as they as
%   simply added to the global wheel stiffness matrix separately.
    kR = cell(ns,1);
    Kcont = cell(ns,1);
%**Compute the absolute angles between the spokes
    thetas = zeros(1,ns);
    t0 = xC(ns+1);
    for j = 1:ns
        thetas(j) = sum(xC(ns+1:ns+j)) - t0;
    end
    
%**Generate the matrices
if usepar
    parfor j = 1:ns
%        fprintf('Building Contact Matrix %d/%d...\n',j,ns); 
        % Rotate element stiffnesses
            [Rmat] = RblockGen(thetas(j),-pi/2);
            kR{j} = Rmat*k0*Rmat';
        % Compute the contact stiffness matrix for the wheel layer
            Xcont = xC(j)*cos(thetas(j)-pi/2);
            Ycont = xC(j)*sin(thetas(j)-pi/2);
            Kcont{j} = contactStiffnessNEW([Xcont;Ycont], sig, Kcont0, ndof, lx, ly, Ex, Ey);
%        clc;
    end
else
    for j = 1:ns
%        fprintf('Building Contact Matrix %d/%d...\n',j,ns); 
    % Rotate element stiffnesses
        [Rmat] = RblockGen(thetas(j),-pi/2);
        kR{j} = Rmat*k0*Rmat';
    % Compute the contact stiffness matrix for the wheel layer
        Xcont = xC(j)*cos(thetas(j)-pi/2);
        Ycont = xC(j)*sin(thetas(j)-pi/2);
        Kcont{j} = contactStiffnessNEW([Xcont;Ycont], sig, Kcont0, ndof, lx, ly, Ex, Ey);
%        clc;
    end
end
        
%--------------------------------------------------------------------------
%   BEGIN COMPUTATION
%--------------------------------------------------------------------------         
%**Initialize the design variable vector        
    rho = HeavisideInv(Vmax,H_beta)*ones(ne, 1);
    
%**Element connectivity filtering
%   NOTE: 1 corresponds to the BARROW layer; 2 is the WHEEL layer
    rho1 = W_filt*rho(1:lne);       
    rho2 = W_filt*rho(lne+1:ne);
    rho_filt = [rho1;rho2];

%**Perform heaviside filtering
    [rho_fH,drHdr] = Heaviside(rho_filt,H_beta);
    
% *************************************************************************
%
% *************************************************************************
close; clc; tic;

C = zeros(ns, 1);
dCdx = zeros(ne, ns);

drHdrB = drHdr(1:lne);
drHdrW = drHdr(lne+1:ne);

if usepar
    parfor i = 1:ns
        % fprintf('Computing compliance case %d/%d...\n',i,ns);
        % Evaluate compliance and compliance sensitivities
        [C(i), dCdX] = ...
            compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                                axle, Ex, Ey, nodeMap, k0, kR{i}, Kcont{i});
        dCdx(:, i) = [((dCdX(1:lne).*drHdrB)'*W_filt)';((dCdX(lne+1:ne).*drHdrW)'*W_filt)'];
        % clc;
    end
else
    for i = 1:ns
    % fprintf('Computing compliance case %d/%d...\n',i,ns);
    % Evaluate compliance and compliance sensitivities
    [C(i), dCdX] = ...
        compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                            axle, Ex, Ey, nodeMap, k0, kR{i}, Kcont{i});
    dCdx(:, i) = [((dCdX(1:lne).*drHdrB)'*W_filt)';((dCdX(lne+1:ne).*drHdrW)'*W_filt)'];
    % clc;
    end
end
    

maxC = max(C);
Cmax = 1.5*maxC; % Maximum allowable compliance

L = ones(lne, 1);
Vol = sum(rho_fH);
dVol(1:lne) = ((L.*drHdrB)'*W_filt)';
dVol(lne+1:ne) = ((L.*drHdrW)'*W_filt)';

%%

%--------------------------------------------------------------------------
%   Initialize parameters used in MMA problem
%--------------------------------------------------------------------------
%**Initial values
    xval = rho;
    f0val = (Vol);
    df0dx = dVol';
    fval = C - Cmax;
    dfdx = dCdx;
%**MMA Parameters
    m = ns;
    epsimin = 1e-7;
    xold1   = xval;
    xold2   = xval;
    xmin    = zeros(ne, 1);
    xmax    = ones(ne, 1);
    low     = xmin;
    upp     = xmax;
    c_mma   = 1000*ones(m, 1);
    d_mma   = 1*ones(m, 1);
    a0_mma  = 1;
    a_mma   = zeros(m, 1);
    outeriter = 0;
    kktnorm = kkttol + 1;
%**Create history arrays
    x_hist = zeros(length(xval),maxoutit);
    C_hist = zeros(maxoutit,1);
    V_hist = zeros(maxoutit,1);
%**Get ready to start optimizing! Print first values
    fprintf('Iter.  KKT norm     max(C)     V/Vmax  %c   %ct     %ct\n',951,916,8747);
    fprintf('----- ---------- ------------ -------- - ------ ------\n');
    T = toc;
    fprintf('%5d %10.4f %12.6f %8.3f %1d %6.2f %6.2f\n',...
        outeriter,kktnorm,maxC,Vol/ne,p_SIMP,T,T/60);
    tic;

%%
while outeriter < maxoutit && kktnorm > kkttol

    outeriter = outeriter+1;
    
    %**Continuation of SIMP parameter (to increasingly penalize elements
    %   with intermediate volume fractions) and of the Heaviside filter
    %   parameter beta (as to create a crisper mesh, with sharper
    %   transitions between void and full material).
    p_SIMP = continuation(p_SIMP,outeriter,1,50,maxSIMP);
    % H_beta = continuation(H_beta,outeriter,50,200,700);
        
    %%%% The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
        mmasub(m,ne,outeriter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval,dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma,H_beta);
    
    %%%% Some vectors are updated:
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;

    % layerwise filtered densities
    rho1 = W_filt*xval(1:lne); 
    rho2 = W_filt*xval(lne+1:ne);
    [rho_fH,drHdr] = Heaviside([rho1;rho2],H_beta);
       
    % SENSITIVITY
    drHdrB = drHdr(1:lne);
    drHdrW = drHdr(lne+1:ne);

    if usepar
        parfor i = 1:ns
            % Evaluate compliance and compliance sensitivities
            [C(i), dCdX] = ...
                compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                                    axle, Ex, Ey, nodeMap, k0, kR{i}, Kcont{i});
            dCdx(:, i) = [((dCdX(1:lne).*drHdrB)'*W_filt)';((dCdX(lne+1:ne).*drHdrW)'*W_filt)'];
        end
    else
       for i = 1:ns
        % Evaluate compliance and compliance sensitivities
        [C(i), dCdX] = ...
            compliance2D_multi(rho_fH, p_SIMP, F_ext, constraints, ...
                                axle, Ex, Ey, nodeMap, k0, kR{i}, Kcont{i});
        dCdx(:, i) = [((dCdX(1:lne).*drHdrB)'*W_filt)';((dCdX(lne+1:ne).*drHdrW)'*W_filt)'];
       end
    end

   
    Vol = sum(rho_fH);   
    dVol(1:lne) = ((L.*drHdrB)'*W_filt)';
    dVol(lne+1:ne) = ((L.*drHdrW)'*W_filt)';
    
    % Objective and constraint functions, sensitivities updated
    f0val = Vol;
    df0dx = dVol';
    fval = C - Cmax;
    dfdx = dCdx;
    
    %%%% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = kktcheck(m,ne,xmma,ymma,zmma,...
        lam,xsi,eta,mu,zet,s,xmin,xmax,df0dx,fval,dfdx',...
        a0_mma,a_mma,c_mma,d_mma);
    
    T = T + toc;
    fprintf('%5d %10.4f %12.6f %8.3f %1d %6.2f %6.2f\n',...
        outeriter,kktnorm,max(C),Vol/ne,p_SIMP,toc,T/60);
    tic;

    %**PLOT COMPONENT TOPOLOGIES
    figure(2);  % BARROW/Carriage layer plotting
        structplot(rho_fH(1:lne),Ex,Ey,lx,ly);
        title(['Barrow | Iteration: ',num2str(outeriter)]);
        animation(1,outeriter) = getframe(2);
    figure(3);  % WHEEL layer plotting
        structplot(rho_fH(lne+1:end),Ex,Ey,lx,ly); 
        title(['Wheel | Iteration: ',num2str(outeriter)]);
        animation(2,outeriter) = getframe(3);
        
    if (outeriter == 100) || (outeriter == 200) || (outeriter == 400)
        outname = ['TO_iter',num2str(outeriter),'_'];
        makepng(animation(1,outeriter),[outname,'Barrow']);
        makepng(animation(2,outeriter),[outname,'Wheel']);
    end
        
        
    %**Store variable and results history
        x_hist(:, outeriter) = xval;
        C_hist(outeriter) = max(C);
        V_hist(outeriter) = Vol;
        
    % Repeat...
end


makepng(animation(1,end),'TO_finalbarrow');
makepng(animation(2,end),'TO_finalwheel');

makegif(animation(1,1:end),'TO_GIFbarrow');
makegif(animation(2,1:end),'TO_GIFwheel');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -----  |   |  \   |  /----  _____  -----  /---\  \   |  /----   %
%   |      |   |  |\  |  |        |      |    |   |  |\  |  |       %
%   ----   |   |  | \ |  |        |      |    |   |  | \ |  \---\   %
%   |      |   |  |  \|  |        |      |    |   |  |  \|      |   %
%   |      \___/  |   \  \____    |    -----  \___/  |   \  ____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [C, dC] = ...
    compliance2D_multi(rho, p, F, cons, axle, Ex, Ey, nodeMap, k0, kR, Kcont)

% function calculates the displacements and compliance of the topology 
% defined by the relative density matrix rho, subject to the applied force F 
% the vector 'constraints' contains the numbers of the constrained 
% degrees of freedom and
% dimx and dimy are the number of elements along the x and y dimensions 
% of the domain
    
    ndof = (Ey+1)*(Ex+1)*2;
    el_dim = 2;

    lne = Ex*Ey; % number of elements per layer
    ne = 2*lne;  % total number of elements

    [K1, K2] = global_stiffness2D_multi(rho, p, kR, k0, el_dim, Ex, Ey, nodeMap);
    
    % define constraints/supports
    ndofGlob = 2*ndof-2; % # of dofs in global system with axle constraint
    freedofs = setdiff(1:ndofGlob, cons);


    Uglob = zeros(ndofGlob, 1);
    Fglob = sparse(ndofGlob, 1);

    % Assemble global stiffness matrix from layer stiffness matrices
    L1dofs = 1:ndof; % layer 1 dofs
    L2dofs = [ndof+1:ndof+axle(1)-1, axle, ndof+axle(1):ndofGlob]; % all layer 2 dofs expressed in global indexing
    Kglob = sparse(ndofGlob, ndofGlob);
    Kglob(L1dofs, L1dofs) = Kglob(L1dofs, L1dofs) + K1;
    Kglob(L2dofs, L2dofs) = Kglob(L2dofs, L2dofs) + K2 + Kcont; % Add stiffness due to contact (friction) with the ground
    
    % Solve linear systems and re-assemble global matrices/vectors
    Fglob(L1dofs) = F;
    Uglob(freedofs) = Kglob(freedofs, freedofs)\Fglob(freedofs, :);
    C = Fglob'*Uglob;
    
    % Extract layer-specific displacement vectors
    U1 = Uglob(L1dofs);
    U2 = Uglob(L2dofs);

    dC = zeros(ne, 1);
    edof = zeros(2*el_dim^2, 1);

    % Compute compliance sensitivities
    for i = 1:lne
        % global degrees of freedom affected by element (i, j) 
        for m = 1:el_dim^2
            edof(2*m - 1) = 2*nodeMap(i,m) - 1;
            edof(2*m) = 2*nodeMap(i,m);
        end
        dC(i,1) = -p*rho(i)^(p-1)*U1(edof).'*k0*U1(edof);
        dC(i+lne,1) = -p*rho(i+lne)^(p-1)*U2(edof).'*kR*U2(edof);
    end
end



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


function Kcont = contactStiffnessNEW(pin, sig, Kcont0, ndof, lx, ly, Ex, Ey)
% Function computes the stiffness contributions due to the wheel's contact
% with the ground.

    Kcont = zeros(ndof,1);   % initialize modified global (layer) stiffness

    shoplot = 0;
    
    dx = lx/Ex;
    dy = ly/Ey;

    kl = Kcont0*ones(2,1);         % link stiffness
    coeff = sqrt(2*pi)*sig;     % coefficient used in Gaussian function evaluation

    if shoplot == 1
        vis = zeros(Ey+1, Ex+1); % Array used for visializing placement of contact stiffness
    end
    
    for j = 1:Ey+1
        
        Y = -ly/2 + (j-1)*dy;       % Y-location of current node
        
        for k = 1:Ex+1
         
            X = -lx/2 + (k-1)*dx;   % X-location of current node
            dist2 = (X - pin(1))^2 + (Y - pin(2))^2; % squared distance between FEM node and spoke tip
            Phi = coeff*exp(-dist2/(2*sig^2));

            if shoplot == 1
                vis(j, k) = vis(j, k) + Phi;
            end

            k_i = Phi*kl;
            
            % Find where the stiffness is to be added to the matrix
            nodeNum = (j-1)*(Ex+1) + k;
            edofs = [2*nodeNum-1, 2*nodeNum];
            % Add it in!
            Kcont(edofs,1) = Kcont(edofs,1) + k_i;
        end
    end
    
    Kcont = spdiags(Kcont,0,ndof,ndof);

    % Plot contact map
    if shoplot == 1
        figure(6); cla; axis equal;
        [Xgrid, Ygrid] = meshgrid(-lx/2:dx:lx/2, -ly/2:dy:ly/2);
        surf(Xgrid, Ygrid, vis); hold on
        plot(pin(1), pin(2), 'm*'); hold on
        view([0 90]);
    end
end


function [nodeMap, X, Y] = createMesh_multi(Ex, Ey, lx, ly)

    % Assumes linear elements

    nel = Ex*Ey; % Number of elements per layer
    nodeMap = zeros(2*nel, 4);

    dx = lx/Ex;
    dy = ly/Ey;

    % Assume both layers are the same with identical numbering

    for j = 1:Ey
        for i = 1:Ex

            ind = (j-1)*Ex + i;

            nodeMap(ind, 3) = (j-1)*(Ex+1) + i;
            nodeMap(ind, 4) = (j-1)*(Ex+1) + i+1;
            nodeMap(ind, 1) = (j)*(Ex+1) + i;
            nodeMap(ind, 2) = (j)*(Ex+1) + i+1;

        end
    end

    for j = 1:Ey+1
        for i = 1:Ex+1

            ind = (j-1)*(Ex+1) + i;
            X(ind, 1) = (i-1)*dx;
            Y(ind, 1) = ly - (j-1)*dy;

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


% Function generates the stiffness matrix for an arbitrary quadrilateral 
% element given the Young's modulus E, Poisson's ratio v, and vectors
% containing the X and Y coordinates of the nodes of the element
% The function uses an isoparametric formulation
function [ke, vol] = Ke2Diso(E, v, X, Y)

    % Solve for stiffness Matrix

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

%%%%%%%%%%%%%%%%%%%%%% Shape Functions %%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates the DERIVATIVES of the shape function
% given the natural coordinates, s, t
function dN = shapeFunctions(s, t)
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

%%%%%%%%%%%%%%%%%%%%%% Strain-Displacement Matrix %%%%%%%%%%%%%%%%%%%%%%%%%
function [B, detJ] = StrainDisp(X_el, Y_el, s, t)

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

%%%%%%%%%%%%%%%%%%%%%% Constitutive Matrix %%%%%%%%%%%%%%%%%%%%%%%%%
function [D, M] = constitutive(E, v)

  % Generate constitutive matrix D given Young's Modulus E, and
  % Poisson's ratio v
  
  D(1, 1) = E/(1-v^2);
  D(1, 2) = v*E/(1-v^2);
  D(2, 1) = D(1, 2);
  D(2, 2) = D(1, 1);
  D(3, 3) = E/(2*(1+v));
  
  M(1, 1) = 1;
  M(1, 2) = -1/2;
  M(2, 1) = M(1, 2);
  M(2, 2) = M(1, 1);
  M(3, 3) = 3;
  
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
%---------------------------------------------------------------------
end


%-------------------------------------------------------
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,H_beta)
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
    epsimin = 10^(-7);
    raa0 = 0.00001;
    albefa = 0.95;
    asyinit = 0.5/(1+H_beta);
    % asyincr = 1.2;
    % asydecr = 0.7;
    asyincr = 1.1;      % Modified by PK
    asydecr = 0.9;      % "            "
    eeen = ones(n,1);
    eeem = ones(m,1);
%    zeron = zeros(n,1);

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
     %  lowmax = xval - 0.01*(xmax-xmin);
     %  uppmin = xval + 0.01*(xmax-xmin);
      lowmax = xval - min(0.01,1/H_beta)*(xmax-xmin);     % Modified by PK
      uppmin = xval + min(0.01,1/H_beta)*(xmax-xmin);     % "
      % uppmax = xval + 10*(xmax-xmin);
      uppmax = xval + 0.05*(xmax-xmin);                      % "
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
%    p0 = zeron;
%    q0 = zeron;
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
%    P = sparse(m,n);
%    Q = sparse(m,n);
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
%-------------------------------------------------------------
end


function structplot(rho,nex,ney,Lx,Ly) %,X,Y)

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


function [xCE] = ellipsegen(ns)

    % Ellipse parameter as computed from optimal shape
    a = 1.5654;
    b = 1.0825;
    c = sqrt(a^2 - b^2);
    e = c/a;
    
    % Angles for equal spoke spacing
    dQ = 2*pi/ns;
    E = 0:dQ:2*pi-dQ;
    f = 2*atan2(tan(E/2),sqrt((1+e)/(1-e)));
    for j = 1:ns
        if f(j) < 0
            f(j) = 2*pi + f(j);
        end
    end

    % Radii and angular coordinates to be sent out
    RE = a.*(1+e.*cos(E));
    QE = [0 diff(f)];
    xCE = [RE,QE]';

end
