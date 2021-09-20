% TOPOLOGY OPTIMIZATION MAIN PROGRAM
%  for 2D planar, 4-node isoparametric quadrilateral elements.
%   Objective function type: COMPLIANCE MINIMIZATION,
%   i.e. f = F'Q
%  Uses Heaviside filtering (as in Guest, 2011) and the MMA as
%   the optimizer, with beta-coninutation eliminated.
%  Assumes linear continuum mechanics FEA equations.
%
% PATRICK KELLEY - UIUC - 2020
%
clear; clc; close all;
%
%**************************************************************
%** USER INPUTS ***********************************************
%**************************************************************

  %** Physical length of the mesh [m], in the X & Y directions,
  %    along with the (constant) thickness in the Z direction
    L_x = 4;
    L_y = 2;
    t = 1;
    
  %** Number of elements in the X and Y directions for the mesh
    nex = 200;
    ney = 100;
    
  %** Forces  
    P = -2;
    % (x,y) location of the applied force
    F_X = L_x;
    F_Y = L_y/2;
    
  %** CONSTITUTIVE RELATIONS ********************************** 
  %* stype : Stress-strain situation: give string 'stress'    *
  %           for plane-stress problems; 'strain' similarly   *
  %**** YM : Young's Modulus, [N/m] (or equivalent)           *
  %**** nu : Poisson's Ratio (usually 0.3 for metals)         *
  %************************************************************
    stype = 'stress';   
    YM = 1;
    nu = 0.3;
    
  %** TOPOLOGY OPTIMIZATION PARAMETERS ************************
  %* rho_min : minimum allowable element density, should be   *
  %              close to zero. Used to prevent singularness  *
  %              of the stiffness matrix K.                   *
  % filt_rad : radius around an element which to consider as  *
  %              the "neighborhood" for filtering. Give in [m]*
  %** H_beta : smoothing parameter for Heaviside Function     *
  %              1 - exp(-x*H_beta) + x*exp(-H_beta). Also    *
  %              controls the asymptotes in the MMA problem   *
  %* max_eta : maximum exponent for SIMP penalty procedure    *
  %              which starts at 1 and, in a continuation     *
  %              method, increases up to max_eta when local   *
  %              convergence (or the maximum number of        *
  %              iterations, see below) is achieved.          *
  %* cont_it : number of iterations (as an upper bound,       *
  %              assuming no local convergence) before the    *
  %              SIMP continuation parameter for SIMP is      *
  %              increased by 1.                              *
  % cont_con : change in objective function value at which    *
  %              (local) convergence is assumed for a given   *
  %              continuation step, allowing the SIMP         *
  %              exponent to be increased by one.             *
  % maxoutit : maximum number of iterations for the optimizer.*
  %************************************************************
    rho_min = 1e-6;
    filt_rad = 0.05;
    H_beta = 500;
    max_eta = 5;
    cont_it = 50;
    cont_con = 1e-3;
    maxoutit = 500; %(max_eta+1)*cont_it;
    
  %** CONSTRAINTS *********************************************
  % vol_frac : maximum allowable volume fraction.             *
  %************************************************************
    vol_frac = 0.5;

    
%--------------------------------------------------------------    
%-- A FEW QUICK, UP-FRONT CALCULATIONS ------------------------
%--------------------------------------------------------------

  % Compute the spatial element scale in the two directions
    dx = L_x/nex;
    dy = L_y/ney;
  % Compute the volume of the element
    ve = dx*dy*t;
  % Normalized total volume
    V0 = nex*ney;
  % Generate constitutive matrix
    D = matmat(YM,nu,stype);
  % Initial SIMP penalty exponent starts off at
    eta = 2;
    eta0 = eta;
    
%-- GRAB MESH
    [nodes,numno,elems,numel,BC,nfree,nfix] = ...
        qmeshr([L_x,L_y],[nex,ney],{'left'});    
    
%-- CREATE LOADS
    loads = sparse(2,numno);
    nP = fnodexy(F_X,F_Y,[L_x,L_y],[nex,ney]);
    loads(2,nP) = P;
    F = reshape(loads,2*numno,1);
    
%-- DEFAULT ELEMENT STIFFNESS MATRIX
    K0 = buildqsm(1,nodes,elems,D,t);

%-- GENERATE THE FILTER MATRIX
    fprintf('Building Filter Matrix...\n');
    W = denfiltmat(filt_rad,elems,nodes);
%     W = eye(numel); % #nofilter
    clc;
    
%-- CREATE INITIAL DESIGN VARIABLES
    x0 = HeavisideInv(vol_frac-1e-3,H_beta)*ones(numel,1);
    

%% ------------------------------------------------------------
% INITIAL COMPUTATIONS
    tic;

  % Perform distance filtering
    x = W*x0;
  % Perform Heaviside filtering -- generate physical variables
    [X,dXdx] = Heaviside(x,H_beta);
  % Compute the compliance and its sensitivity
    [C,dCdX] = complisub(X,elems,K0,F,nfree,eta,rho_min);
  % Compute the volume constraint
      % Compute the current volume fraction for the constraint
      %  g(x) = sum(x)/V0 - vol_frac <= 0
        V = sum(X)/V0 - vol_frac;
      % The sensitivities of the volume are linear in variables
        dVdX = ones(size(X))/V0;
  % Defilter the sensitivities
    dCdx = W'*((dCdX).*(dXdx));
    dVdx = W'*((dVdX).*(dXdx));
    
    T = toc;
    
    figure(1);
    structplot(X,nex,ney,L_x,L_y); axis off;
    
%% ------------------------------------------------------------
% Initialize parameters used in MMA problem
    xval = x0;                      % Initial variable values
    m = 1;                          % No. of constraints
    xold1   = xval;                 % Previous variables are the initial
    xold2   = xval;                 % "                                "
    xmin    = 0*ones(numel, 1);     % Minimum allowable values
    xmax    = 1.0*ones(numel, 1);  % Maximum "              "
    low     = xmin;
    upp     = xmax;
    c_mma   = 1000*ones(m, 1);
    d_mma   = 1*ones(m, 1);
    a0_mma  = 1;
    a_mma   = 0*ones(m, 1);
    outeriter = 0;                  % Initial iteration number
%     maxoutit  = 300;                % Total number of iterations to run
    kkttol  = 2.0e-4;               % KKT condition tolerance
    kktnorm = kkttol + 1;           % Initial KKT gradient norm
    f0val = C;                      % Initial work of the shape
    df0dx = dCdx;                   % Initial gradient
    fval = V;                       % Constraint value
    dfdx = dVdx;                    % Constraint gradient
    
    
%% ------------------------------------------------------------
% Print the table header for the iterations
    fprintf('Iter.  KKT norm      C(x)       V/V0   %c   %ct     %ct\n',951,916,8747);
    fprintf('----- ---------- ------------ -------- - ------ ------\n');
    fprintf('%5d %10.4f %12.6f %8.4f %1d %6.3f %6.2f\n',outeriter,kktnorm,C,sum(X)/V0,eta,T,T/60);
    tic;
    
%% ------------------------------------------------------------
%-- OPTIMIZE! -------------------------------------------------
    
cont_flag = 0;

while outeriter < maxoutit && kktnorm > kkttol
    % Increase the iteration number
        outeriter = outeriter + 1;    
        
    if outeriter < 3 || (cont_flag > 0 && cont_flag < 3)
        asyflag = 1;
        cont_flag = cont_flag + 1;
    else
        asyflag = 0;
        cont_flag = 0;
    end
        
    % The MMA subproblem is solved at the point xval:
    [xmma,ymma,zmma,lam,xsi,etam,mu,zet,s,low,upp] = ...
        mmasub(m,numel,outeriter,xval,xmin,xmax,xold1,xold2,...
        f0val,df0dx,fval',dfdx',low,upp,a0_mma,a_mma,c_mma,d_mma,H_beta*asyflag);
    
    % Some vectors are updated:
        xold2 = xold1;
        xold1 = xval;
        xval  = xmma;
        
    % Evaluate objective function, constraints, sensitivities
      % Perform distance filtering
        x = W*xval;
      % Perform Heaviside filtering -- generate physical variables
        [X,dXdx] = Heaviside(x,H_beta);
      % Compute the compliance and its sensitivity
        [C,dCdX] = complisub(X,elems,K0,F,nfree,eta,rho_min);
      % Compute the volume constraint
        V = sum(X)/V0 - vol_frac;
      % Defilter the sensitivities
        dCdx = W'*((dCdX).*(dXdx));
        dVdx = W'*((dVdX).*(dXdx));
        
    % Update the values
        f0val = C;                      % Update the objective function value
        df0dx = dCdx;                   % Update the variable gradient      
        fval = V;                       % Check the perimeter constraint value
        dfdx = dVdx;                    % Update perimeter constraint gradient
        
    % The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,numel,xmma,ymma,zmma,lam,xsi,etam,mu,zet,...
        s,xmin,xmax,df0dx,fval',dfdx',a0_mma,a_mma,c_mma,d_mma);
    
    % Continuation: gradual increase the SIMP exponent
        if ((kktnorm < cont_con) || (outeriter > (eta-eta0+1)*cont_it)) && (eta ~= max_eta)
            eta = eta+1;
            cont_flag = 1;
        end
    
    % Output Information
        T = T + toc;                    % Keep track of the total time
        fprintf('%5d %10.4f %12.6f %8.4f %1d %6.3f %6.2f\n',...
            outeriter,kktnorm,C,sum(X)/V0,eta,toc,T/60);
        tic;                            % Reset the timer
        
    % Plot Design
        figure(1);
        structplot(X,nex,ney,L_x,L_y); axis off;
        title(['Iteration: ',num2str(outeriter)]);
        animation(outeriter) = getframe();
        
    % Save the current status to the history arrays    

    % Repeat...
end    

%-- Inform the user of the optimizer's reason for stopping.
if kktnorm < kkttol
    fprintf('\nKKT conditions satisfied to within tolerance!\n');
else
    fprintf('\nMaximum optimizer iterations reached.\n');
end
    
    
structplot(X,nex,ney,L_x,L_y,'CantBeamFiltHeavi.png'); axis off;
title('Final Result');
        
    




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   -----  |   |  \   |  /----  _____  -----  /---\  \   |  /----   %
%   |      |   |  |\  |  |        |      |    |   |  |\  |  |       %
%   ----   |   |  | \ |  |        |      |    |   |  | \ |  \---\   %
%   |      |   |  |  \|  |        |      |    |   |  |  \|      |   %
%   |      \___/  |   \  \____    |    -----  \___/  |   \  ____/   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [E] = matmat(E0,nu,stype)
%MATMAT   Generate material matrix for an element.
%   
%   INPUTS:
%    E0......Young's Modulus
%    nu......Poisson's Ratio
%    stype...String of the problem type, either 2-D plane-stress,
%             2-D plane-strain, or general 3-D.
%             Input examples: 'pss', 'stress', 'plane stress'
%                             'psn', 'strain', 'plane strain'
%             Do not supply this input if 3-D.
%   OUTPUT:
%    E.......Material Matrix, either 3x3 for 2-D or 6x6 for 3-D
%
% PATRICK KELLEY - UIUC - 2020

    if nargin < 3
        stype = '';
    end

    switch stype
        case {'pss', 'stress', 'plane stress'}
            E = (E0/(1-nu^2)).* ...
                [1 nu 0;
                 nu 1 0;
                0 0 (1-nu)/2];
        case {'psn', 'strain', 'plane strain'}
            % Plane-strain case
            E = (E0/((1+nu)*(1-2*nu))).* ...
                [1-nu nu 0;
                 nu 1-nu 0;
                0 0 (1-nu)/2];
        otherwise
            % General 3-D case
            E = (E0/((1+nu)*(1-2*nu))).* ...
                [1-nu nu nu 0 0 0;
                 nu 1-nu nu 0 0 0;
                 nu nu 1-nu 0 0 0;
                 0 0 0 0.5-nu 0 0;
                 0 0 0 0 0.5-nu 0;
                 0 0 0 0 0 0.5-nu];
    end
end


function [nodes,numno,elems,numel,BC,nfree,nfix,loads,PF] ...
    = qmeshr(rdim,nn,fixno)
%QMESHR   Create a quadrilateral mesh in a rectangular
%          domain.
%
%  INPUTS:
%   rdim....Vector of the rectangular domain's dimensions [m],
%            a 1x2 (or 2x1) vector of the form [L_x L_y]. 
%   nn......Number of elements across the meshes dimensions,
%            a 1x2 (or 2x1) vector of the form [nnx nny] with
%            nnx being the no. of elements in the x-direction.
%   fixno...4x(-) array giving horizontal or vertical lines
%            along which fixed conditions are specified. A 
%            column represents [x0 y0 xf yf]'.
%
%  OUTPUTS:
%   nodes...2x({nnx+1}*{nny+1})
%
%   numno...Total number of nodes ( = {nnx+1}*{nny+1}).
%   elems...4x(nnx*nny) array, where the columns are the 
%            nodal numbers of that element as labeled 
%            in a CCW scheme.
%   numel...Total number of elements ( = nnx*nny).
%   BC......2x({nnx+1}*{nny+1}) array with zeros indicating
%            that dof is free, and a 1 when fixed.
%   nfree...Vector of free DOF.
%   nfix....Vector of fixed DOF.
%
% PATRICK KELLEY - UIUC - 2020

%=============================================================
% NODE LOCATIONS AND CONNECTIVITY
    [nodes,elems] = rectmesh2d(rdim,nn);
    numno = length(nodes);      % Total number of nodes
    numel = size(elems,2);      % Total # of elements
    nsd = 2;
    
%=============================================================
% BOUNDARY CONDITIONS
    
    % 
    if iscell(fixno) == 1
        nnx = nn(1);    nny = nn(2);
        bx = [];
        by = [];
        for j = 1:size(fixno,2)
           switch fixno{j}
           % Cases where both are fixed
               case {'top','Top','TOP','t','T'}
                   bx = [bx ((nnx+1)*(nny)+1):((nnx+1)*(nny+1))];
                   by = bx;
               case {'bottom','Bottom','BOTTOM','b','B'}
                   bx = [bx 1:(nnx+1)];
                   by = bx;
               case {'left','Left','LEFT','l','L'}
                   bx = [bx 1:(nnx+1):(((nnx+1)*nny)+1)];
                   by = bx;
               case {'right','Right','RIGHT','r','R'}
                   bx = [bx (nnx+1):(nnx+1):((nnx+1)*(nny+1))];
                   by = bx;
           % Cases where only x is fixed
               case {'xtop','xTop','xTOP','xt','xT'}
                   bx = [bx ((nnx+1)*(nny)+1):((nnx+1)*(nny+1))];
               case {'xbottom','xBottom','xBOTTOM','xb','xB'}
                   bx = [bx 1:(nnx+1)];
               case {'xleft','xLeft','xLEFT','xl','xL'}
                   bx = [bx 1:(nnx+1):(((nnx+1)*nny)+1)];
               case {'xright','xRight','xRIGHT','xr','xR'}
                   bx = [bx (nnx+1):(nnx+1):((nnx+1)*(nny+1))];
           % Cases where only x is fixed
               case {'ytop','yTop','yTOP','yt','yT'}
                   by = [by ((nnx+1)*(nny)+1):((nnx+1)*(nny+1))];
               case {'ybottom','yBottom','yBOTTOM','yb','yB'}
                   by = [by 1:(nnx+1)];
               case {'yleft','yLeft','yLEFT','yl','yL'}
                   by = [by 1:(nnx+1):(((nnx+1)*nny)+1)];
               case {'yright','yRight','yRIGHT','yr','yR'}
                   by = [by (nnx+1):(nnx+1):((nnx+1)*(nny+1))];
               otherwise
                   fprintf('Not a possible fix side\n');
                   return;
           end
            
        end
        
        
    else
        nBCline = size(fixno,2);    % No. of lines is cols of fixno

        for j = 1:nBCline

        end
    end
    
    BC = fastBC(bx,by,nodes);

    BCv = reshape(BC,1,length(BC)*nsd);
    fdof = length(BCv)-sum(BCv);  % Number of free degrees of freedom

    % Set up which DOFs are of interest and which are fixed
    nfree = zeros(1,fdof);
    nfix = zeros(1,length(BCv)-fdof);
    m = 1; ki = 1;
    for j = 1:length(BCv)
        if BCv(j) == 1
            nfix(m) = j;
            m = m+1;
        else
            nfree(ki) = j;
            ki = ki+1;
        end
    end
end


function [nodes,elems] = rectmesh2d(dim,es)
%RECTMESH2D   Quickly create a mesh ovre a rectangular region.
%
%  INPUTS:
%   dim.....2x1 vector with the first element giving length
%            and the second width
%   es......2x1 vector specifying number of elements in the
%            x and y directions respectively
%  OUTPUTS: 
%   nodes...node location array, 2x(numno)
%   elems...element connectivity array, 4x(numel)
%
% PATRICK KELLEY - UIUC - 2020

    L = dim(1);
    W = dim(2);
    ex = es(1);
    ey = es(2);

    xspc = L/ex;
    yspc = W/ey;
    nelem = ex*ey;
    
    nodes = zeros(2,(ex+1)*(ey+1));
    elems = zeros(4,nelem);
    
    % Store the locations of all the nodes
    k = 1;
    for j = 0:ey
        for i = 0:ex
            nodes(1,k) = xspc*i;
            nodes(2,k) = yspc*j;
            k = k+1;
        end
    end

    k = 0;
    % Create the element interconnecteness table
    for j = 1:nelem
        
        if (mod(j-1,ex) == 0) && j ~= 1
            k = k+1;
        end
        
        elems(1,j) = j+k;
        elems(2,j) = j+1+k;
        elems(3,j) = j+ex+2+k;
        elems(4,j) = j+ex+1+k;
    end
end


function BC = fastBC(fixx,fixy,nodes)
%FASTBC   Quickly create the boundary condition array,
%          knowing which nodes are fixed.
%
%   INPUTS:
%    fixx....Vector of node numbers fixed in the X-direction
%    fixy...."                                 " Y-direction
%   OUTPUTS:
%    BC......Matrix the same size as node array, with ones
%             present where a fixed BC is present, 0 elsewise
%
% PATRICK KELLEY - UIUC - 2020

    BC = zeros(size(nodes));
    BC(1,fixx) = ones(1,length(fixx));
    BC(2,fixy) = ones(1,length(fixy));

end


function n = fnodexy(x,y,rdim,nn)
%FNODEXY   Find the node number of a specific (x,y) location.
%           Assumes a rectangular mesh with [length, height]
%           rdim and node discretization nn and that the (x,y)
%           location lands directly on a node.
%
%   INPUTS:
%    x......x-coordinate
%    y......y-coordinate
%    rdim...[length, width] of the rectangular mesh
%    nn.....vector with no. of elements in x and y
%   OUTPUT:
%    n......node number at (x,y)
%
% PATRICK KELLEY - UIUC - 2020

    Lx = rdim(1);   Ly = rdim(2);
    nnx = nn(1);    nny = nn(2);
    dx = Lx/nnx;    dy = Ly/nny;


    xx = x/dx + 1;
    yy = y/dy;
    
    n = (nnx+1)*yy + xx;

end


function ke = buildqsm(e,nodes,elems,D,te)
%BUILDQSM   Build a stiffness matrix for a quadrilateral element.
%
%   Assumes a 4-node quadrilateral element with 
%    4-pt, 2x2 Gaussian Quadrature, linear elasticity,
%    and uniform thickness.
%
%   INPUTS:
%    e.......element to evaluate
%    nodes...spatial location of the nodes, 2 x numno
%    elems...element connectivity matrix, 4 x numels
%    D.......element material properties matirx, 3x3
%    te......thickness of the element; default is unity
%   OUTPUTS:
%    ke......8x8 element stiffness matrix
%
% PATRICK KELLEY - UIUC - 2020
    
    if size(nodes,1) < size(nodes,2)
        nodes = nodes';
    end
    
    % Find the locations of the nodes in the current element only --
    %  that's all we need to work with
    nodee = nodes(elems(:,e),:);
    
    % Create the values of the gaussian integration points
    % Note that the point weights are all 1.0
    xi_int = (1/sqrt(3))*[-1 1 1 -1];       %  xi integration pt coords
    etaint = (1/sqrt(3))*[-1 -1 1 1];       % eta "                   "
    
    % Evaluate the shape functions and their derivatives at the integration
    %  points. Values along COLUMNS (same row) of N are at same int. pt.
    %[~,Ndxi,Ndeta] = shapefns(xi_int,etaint,4);     % Shape fns at int pts 
    Ndxi = zeros(4,4);
    Ndeta = zeros(4,4);
    for i = 1:4
        for k = 1:4
            Ndxi(i,k) = dNquad(xi_int(i),etaint(i),k,0);
            Ndeta(i,k) = dNquad(xi_int(i),etaint(i),k,1);
        end
    end


    % Preallocate the element stiffness matrix
    ke = zeros(8,8);
    
    % Find the updated shape fcn derivatives in the current config
    [Ndx,Ndy,jj] = Nder(Ndxi,Ndeta,nodee(:,1),nodee(:,2),4,4);
    
    % Loop over the integration points. Note that the weights for
    %  2x2 quadrature are 1; they are not included.
    for i = 1:4    
        B = buildB(Ndx,Ndy,i,4);    % Create the strain-displacement matrix    
        
        ke = ke + (te*jj).*(B'*D*B); 
    end
 
end


% function [N,Ndxi,Ndeta] = shapefns(xi,eta,numno)
% %SHAPEFNS   Evaluate element shape functions and their
% %   derivatives at integration point vectors xi and eta.
% %
% %   Assumes the following isoparametric mapping:
% %
% %   (4) o-----------o (3)
% %       |           |
% %       |           |
% %       |           |
% %       |           |
% %   (1) o-----------o (2)
% %
% %   Values along COLUMNS (i.e same row) of N are at the same int. pt.   
% %
% % PATRICK KELLEY - UIUC - 2019
%     
%     nint = length(xi);
% 
%     N = zeros(nint,numno);
%     Ndxi = zeros(nint,numno);
%     Ndeta = zeros(nint,numno);
% 
%     for i = 1:nint
%         for k = 1:numno
%             N(i,k) = Nquad(xi(i),eta(i),k);
%             Ndxi(i,k) = dNquad(xi(i),eta(i),k,0);
%             Ndeta(i,k) = dNquad(xi(i),eta(i),k,1);
%         end
%     end
%     
% end


function N = Nquad(xi,eta,node)
%NQUAD   Create the shape function polynomials, at
%         the specific (xi,eta) location
%
% PATRICK KELLEY - UIUC - 2019

    xiv = [-1 1];
    etav = [-1 1];

    arel = [1 1 1;
            2 2 1;
            3 2 2;
            4 1 2];
        
    b = arel(node,2);
    c = arel(node,3);
    
    lb = 1;
    lc = 1;
    
    for i = 1:2
       if (i == b)
           continue;
       end
       lb = lb*(xi - xiv(i))/(xiv(b) - xiv(i));
    end
    
    for i = 1:2 
       if (i == c)
           continue;
       end
       lc = lc*(eta - etav(i))/(etav(c) - etav(i)); 
    end

    N = lb*lc;
    
end


function NP = dNquad(xi,eta,node,dir)
%DNQUAD   Find derivatives of the shape functions at a point
%          using finite difference
%   Supply dir as 0 for N_{node,xi}
%                 1 for N_{node,eta}
%
% PATRICK KELLEY - UIUC - 2019
    
    h = 1e-10;      % Finite difference step
    NP = (Nquad(xi + 0.5*(1-dir)*h, eta + 0.5*dir*h,node) - ...
          Nquad(xi - 0.5*(1-dir)*h, eta - 0.5*dir*h,node))/h;

end


function [Ndx, Ndy,jj] = Nder(Ndxi,Ndeta,x,y,nint,numno)
%NDER   Find the derivatives of the shape functions at the
%        integration points in the current configuration.
%
% PATRICK KELLEY - UIUC - 2019

    Ndx = zeros(nint,numno);   Ndy = zeros(nint,numno);
    
    for i = 1:nint    % Loop over int pts
        xdxi = 0;  ydxi = 0;
        xdeta = 0; ydeta = 0;
        
        for j = 1:numno   % Over nodes            
            xdxi = xdxi + Ndxi(i,j)*x(j);
            xdeta = xdeta + Ndeta(i,j)*x(j);
            ydxi = ydxi + Ndxi(i,j)*y(j);
            ydeta = ydeta + Ndeta(i,j)*y(j);
        end
        
        jj = xdxi*ydeta - xdeta*ydxi;
        
        for j = 1:numno
            Ndx(i,j) = (Ndxi(i,j)*ydeta - Ndeta(i,j)*ydxi)/jj;
            Ndy(i,j) = -(Ndxi(i,j)*xdeta - Ndeta(i,j)*xdxi)/jj;
        end
  
    end
end


function B = buildB(Ndx,Ndy,ipt,numno)
%BUILDB    Generate the shape function derivative matrix.
%           Assumes a 2-D element.
%
%   INPUTS:
%    Ndx,Ndy --> Matrices where each the i,j location is  
%                 N_{node j,dir}(ipt_i), i.e. the derivative of
%                 the jth node's shape function at integration
%                 point i (use function SHAPEFNS and NDER to 
%                 generate these at supplied integration points) 
%    ipt ------> Integration point to evaluate at
%    numn0 ----> Number of nodes in the element
%
%   OUTPUT:
%    B --------> strain-displacement matrix, size 3x(2*numno)
%
% PATRICK KELLEY - UIUC - 2019

    B = zeros(3,numno*2);
    for a = 1:numno
        B(1,2*a-1) = Ndx(ipt,a);
        B(2,2*a) = Ndy(ipt,a);
        B(3,2*a-1) = Ndy(ipt,a);
        B(3,2*a) = Ndx(ipt,a);
    end

end


function W = denfiltmat(rad,elems,nodes)
%DENFILTMAT   Create a density filter matrix for a
%              rectangular, quadrilateral-element mesh.
%
%  INPUTS:
%   rad.....Filtering radius, [m]
%   elems...Element numbering array, 4 x numel
%   nodes...(x,y) location array for all nodes
%
%  OUTPUT:
%   W.......Density filter matrix. All rows sum to 1.
%
% PATRICK KELLEY - UIUC - 2020

    numel = length(elems);
   
    % Compute the centroid locations of all the elements
    X_cent = zeros(numel,1);
    Y_cent = zeros(numel,1);
    
    for i = 1:numel
        X_cent(i) = sum(nodes(1,elems(:,i)'))/4;
        Y_cent(i) = sum(nodes(2,elems(:,i)'))/4;
    end
    
    % Create the filter matrix
    W = zeros(numel,numel);
    
    for i = 1:numel
        weightsum = 0;
        
        for j = 1:numel
            dij = sqrt((X_cent(i) - X_cent(j))^2 + (Y_cent(i) - Y_cent(j))^2);
%             coeff = rad - dij;
            coeff = 1 - dij/rad;
            if coeff > 0
                W(i,j) = coeff;
                weightsum = weightsum + coeff;
            end
        end
        
        W(i,:) = W(i,:)/weightsum;
    end

    W = sparse(W);
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
% Find the inverse of the Heaviside function
    f = @(x) exp((1-x)*beta) - x - exp(beta)*(1-v0);
    Hi = fzero(f,1/beta);
end


function [C,dCdX] = complisub(X,elems,k0,F,nfree,eta,rho_min)
%COMPLISUB    Compliance FEA subproblem
%
%  INPUTS:
%   X.........Physical design variables (element density fraction) vector
%   elems.....4 x numel array with nodal numbers for element j in column j
%   K0........Default stiffness matrix for an element
%   F.........Global force vector
%   nfree.....Array giving numbers of free degrees of freedom
%   eta.......SIMP exponent parameter
%   rho_min...Minimum acceptable element density
%  OUTPUT:
%   C.........Compliance of the structure
%   dCdX......Sensitivity vector of the compliance WRT physical variables
%
% PATRICK KELLEY - UIUC - 2020

    numel = length(elems);
    ndof = length(F);
    Q = zeros(ndof,1);

    % Assemble the global stiffness matrix
    K = Kassem(k0,elems,numel,X,eta,rho_min);
  
    % Solve for the displacements
    Q(nfree) = K(nfree,nfree)\F(nfree);

    % Evaluate the objective function
    C = F'*Q;
    
    % Compute the sensitivities using the direct method
    dCdX = zeros(numel,1);
    kdof = zeros(8,1);
    
    for j = 1:numel
        for m = 1:4
            kdof(2*m-1) = 2*elems(m,j)-1;
            kdof(2*m) = 2*elems(m,j);
        end
        dCdX(j) = (-eta*X(j)^(eta-1))*((Q(kdof)'*k0)*Q(kdof));
    end
end


function K = Kassem(k0,elems,numel,X,eta,rho_min)
%KASSEM
%
%
% PATRICK KELLEY


    k0v = reshape(k0,64,1);
    
    % ASSEMBLE GLOBAL MATRIX
    I = zeros(64,numel);
    J = zeros(64,numel);
    Kvec = zeros(64,numel);
    kdof = zeros(8,1);
    
    if nargin > 3
    % Assembling for topology optimization
        for j = 1:numel
            % Compute the global indices of the current element
            for m = 1:4
                kdof(2*m-1) = 2*elems(m,j)-1;
                kdof(2*m) = 2*elems(m,j);
            end
            % Store the appropriate indexing
            I(:,j) = reshape(repmat(kdof,1,8)',64,1);
            J(:,j) = repmat(kdof,8,1);
            % Store the element stiffness in the abbreviated matrix
            Kvec(:,j) = Kvec(:,j) + k0v*(X(j)^eta + rho_min);
        end
    else
    % Assembling a standard FEA problem, no TO
        for j = 1:numel
            % Compute the global indices of the current element
            for m = 1:4
                kdof(2*m-1) = 2*elems(m,j)-1;
                kdof(2*m) = 2*elems(m,j);
            end
            % Store the appropriate indexing
            I(:,j) = reshape(repmat(kdof,1,8)',64,1);
            J(:,j) = repmat(kdof,8,1);
            % Store the element stiffness in the abbreviated matrix
            Kvec(:,j) = Kvec(:,j) + k0v;
        end
    end
    
    K = sparse(I,J,Kvec);
    
end


%-------------------------------------------------------
%    This is the file mmasub.m
%
function [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = ...
mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2, ...
f0val,df0dx,fval,dfdx,low,upp,a0,a,c,d,Hbeta)
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
%  Hbeta = Smoothing parameter for Heaviside smoothing. (Added by Pat
%           Kelley, 2020, as suggested in Guest-2011.
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
    albefa = 0.1;
    asyinit = 0.5/(1+Hbeta);
    asyincr = 1.2;
    asydecr = 0.7;
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
        lowmax = xval - 0.01*(xmax-xmin);
        uppmin = xval + 0.01*(xmax-xmin);
%         lowmax = xval - min(0.01,1/Hbeta)*(xmax-xmin);     % Modified by PK
%         uppmin = xval + min(0.01,1/Hbeta)*(xmax-xmin);     % "
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
    % p0 = zeron;
    % q0 = zeron;
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
%-------------------------------------------------------------

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


function structplot(rho,nex,ney,Lx,Ly,fname)

    rrho = zeros(ney, nex);

    for j = 1:ney
        for i = 1:nex
            ind = (j-1)*nex + i;
            rrho(ney+1-j, i) = rho(ind);
            % xBlock(j, i) = xVec(ind);
        end
    end
    
    imagesc([0 Lx],[0 Ly],1-rrho);    
    colormap('gray');
    caxis([0 1]);

    if nargin > 5
        imwrite(1-rrho,fname)
    end


    grid off;
    view(0,90);
    xlim([-0.05 1.05]*Lx);
    ylim([-0.05 1.05]*Ly);
    axis equal;
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
