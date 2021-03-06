function [ P U V Xu Yu Xv Yv Xp Yp numYCells numXCells d] = StokesStaggered(g, numYCells, p0, mu, toGraph, height, width, R, L)

% Full Staggered Grid example
%
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-
%   -v-     -v-     -v-     -v-     -v-
%   -p- -u- -p- -u- -p- -u- -p- -u- -p-

% Note that when viewed separately, the p-grid is a 5x5 while the u-grid 
% is 4x5 and the v-grid is 5x4... let d be the spacing of each grid

% For each interior u node, we still use the Stokes equation:
% -(p_xp - p_xm)/d + mu*(u_xpp - 2*u + u_xmm)/d^2 
%       + mu*(u_ypp - 2*u + u_ymm)/d^2 = 0

% For each interior v node, we will eventually use the Stokes equation,
% but to keep things simple, let's use that we know v=0:
% v = 0

% For each interior p node, we will eventually use the full divergence-free
% condition, but to keep things simple, let's use that we know v=0 (implies
% grad v = 0):
% -(u_xp - u_xm)/d = 0

%parameters passed in
% p0 = 0;
% mu = 1.0;
% g = 9.8;
% 
% height = 1;
% width = 10;
% 
% numYCells = 10; 

% Create grids for indexing

d = height/(numYCells-1);

%numXCells = (width / d) + 1;
numXCells = numYCells;

pInd = zeros(numXCells,numYCells);
uInd = zeros(numXCells-1,numYCells);
vInd = zeros(numXCells,numYCells-1);

count = 1;
for i=1:numXCells
    for j=1:numYCells
        pInd(i,j) = count;
        count = count + 1;
    end
end

for i=1:numXCells-1
    for j=1:numYCells
        uInd(i,j) = count;
        count = count + 1;
    end
end

for i=1:numXCells
    for j=1:numYCells-1
        vInd(i,j) = count;
        count = count + 1;
    end
end

[Xp,Yp] = meshgrid(linspace(0,width,numXCells),linspace(-height/2,height/2,numYCells));
[Xu,Yu] = meshgrid(linspace(0.5*d,width-0.5*d,numXCells-1),linspace(-height/2,height/2,numYCells));
[Xv,Yv] = meshgrid(linspace(0,width,numXCells),linspace(-height/2 + .5*d,height/2-0.5*d,numYCells-1));

[X, Y] = meshgrid(linspace(0.5*d,width-0.5*d,numXCells-1), linspace(0.5*d,height-0.5*d,numYCells-1));

[F1u F2u] = force(Xu, Yu, mu, width, p0, g, R, d, L);
[F1v F2v] = force(Xv, Yv, mu, width, p0, g, R, d, L);


count = count-1;  %  Count is 3M^2?
%A = zeros(count); 
A = sparse([],[],[],count,count,5*count);
b = zeros(count,1);

% Iterate through p-grid and set equations
for i=1:numXCells
    for j=1:numYCells 
        
        if (i == 1) % left boundary (use Dirichlet BC)
            
            A(pInd(i,j),pInd(i,j)) = 1.0;
            b(pInd(i,j)) = 0;
            
        elseif (i == numXCells) % right boundary (use Dirichlet BC)
            
            A(pInd(i,j),pInd(i,j)) = 1.0;
            b(pInd(i,j)) = 0;
            
        elseif (j == 1)
            
            A(pInd(i,j),pInd(i,j)) = 1.0;
            b(pInd(i,j)) = 0;
 
        elseif (j == numYCells) 
            
            A(pInd(i,j),pInd(i,j)) = 1.0;
            b(pInd(i,j)) = 0;
            

        else % interior point (use divergence-free condition)
            
            A(pInd(i,j),uInd(i,j)) = 1.0/d;
            A(pInd(i,j),uInd(i-1,j)) = -1.0/d;
            
            A(pInd(i,j),vInd(i,j)) = 1.0/d;
            A(pInd(i,j),vInd(i,j-1)) = -1.0/d;
            
            b(pInd(i,j)) = 0.0;
            
        end
        
    end
end
            
% Iterate through u-grid and set equations
for i=1:numXCells-1
    for j=1:numYCells
               
        if (j == 1) % bottom boundary (use Dirichlet BC)
            
            A(uInd(i,j),uInd(i,j)) = 1.0;
            b(uInd(i,j)) = 0.0;
            
        elseif (j == numYCells) % top boundary (use Dirichlet BC)
            
            A(uInd(i,j),uInd(i,j)) = 1.0;
            b(uInd(i,j)) = 0.0;
            
        elseif (i == 1) % left edge (use modified Stokes)

            %This takes advantage of the fact that u_x = 0 at the boundary
            %ie (u_i - u_i-1) / d = 0
            
            A(uInd(i,j),pInd(i+1,j)) = -1.0/d;
            A(uInd(i,j),pInd(i,j)) = 1.0/d;
            A(uInd(i,j),uInd(i,j)) = -mu*3.0/d^2;
            A(uInd(i,j),uInd(i+1,j)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j+1)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j-1)) = mu*1.0/d^2;
            b(uInd(i,j)) = -F1u(j,i);

            
        elseif (i == numXCells-1) % right edge (use modified Stokes)
            
            A(uInd(i,j),pInd(i+1,j)) = -1.0/d;
            A(uInd(i,j),pInd(i,j)) = 1.0/d;
            A(uInd(i,j),uInd(i,j)) = -mu*3.0/d^2;
            A(uInd(i,j),uInd(i-1,j)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j+1)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j-1)) = mu*1.0/d^2;
            b(uInd(i,j)) = -F1u(j,i);
            
            
        else % interior point (use Stokes)
            
            A(uInd(i,j),pInd(i+1,j)) = -1.0/d;
            A(uInd(i,j),pInd(i,j)) = 1.0/d;
            A(uInd(i,j),uInd(i,j)) = -mu*4.0/d^2;
            A(uInd(i,j),uInd(i+1,j)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i-1,j)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j+1)) = mu*1.0/d^2;
            A(uInd(i,j),uInd(i,j-1)) = mu*1.0/d^2;
            b(uInd(i,j)) = -F1u(j,i);
            
        end
        
    end
end

% Iterate through v-grid and set equations
for i=1:numXCells
    for j=1:numYCells-1
        
        %%% DONE: need to use modified stokes at bottom.
        %%% The "modified" part is that the normal stencil uses v(i,0).
        %%% The value of this can be determined similar to described in
        %%% previous comments
                              
        if (j == 1) % bottom boundary
            
            if (i == 1) %bottom left corner

                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*4.0/d^2;
                A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;
                b(vInd(i,j)) = -F2v(j,i);
            
            elseif (i == numXCells)  % bottom right corner
                
                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*4.0/d^2;
                A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;

                b(vInd(i,j)) = -F2v(j,i);
                
            else  % just the bottom edge
                
                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*5.0/d^2;
                A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;

                b(vInd(i,j)) = -F2v(j,i);                
                
            end

        %%% DONE: need to use modified stokes at top.
        %%% See above comments
        elseif (j == numYCells-1) % top boundary
            
            if(i == 1)  %top left corner
                            
                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*4.0/d^2;
                A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
                b(vInd(i,j)) = -F2v(j,i);
                
            elseif( i == numXCells)  %top right corner
                
                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*4.0/d^2;
                A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
                b(vInd(i,j)) = -F2v(j,i);

            else % just the top edge
                
                A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
                A(vInd(i,j),pInd(i,j)) = 1.0/d;
                A(vInd(i,j),vInd(i,j)) = -mu*5.0/d^2;
                A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
                A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
                b(vInd(i,j)) = -F2v(j,i);
                
            end
            
        elseif (i == 1) % left edge (use modified Stokes)

            A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
            A(vInd(i,j),pInd(i,j)) = 1.0/d;
            A(vInd(i,j),vInd(i,j)) = -mu*3.0/d^2;
            A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
            b(vInd(i,j)) = -F2v(j,i);

            
        elseif (i == numXCells) % right edge (use modified Stokes)
            
            A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
            A(vInd(i,j),pInd(i,j)) = 1.0/d;
            A(vInd(i,j),vInd(i,j)) = -mu*3.0/d^2;
            A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
            b(vInd(i,j)) = -F2v(j,i);
            
            
        else % interior point (use Stokes)
            
            A(vInd(i,j),pInd(i,j+1)) = -1.0/d;
            A(vInd(i,j),pInd(i,j)) = 1.0/d;
            A(vInd(i,j),vInd(i,j)) = -mu*4.0/d^2;
            A(vInd(i,j),vInd(i+1,j)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i-1,j)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j+1)) = mu*1.0/d^2;
            A(vInd(i,j),vInd(i,j-1)) = mu*1.0/d^2;
            b(vInd(i,j)) = -F2v(j,i);
        
        end
    end
end

%sol = pinv(A)*b; % this is needed until we use the vertical velocity to couple everything together
A = sparse(A); % gives Matlab a way to solve this more efficiently
sol = A\b;


% Assign solution vector to grids
% Pstaggered = sol(pInd);
% Ustaggered = sol(uInd);
% Vstaggered = sol(vInd);

% Plots grids

%U = interp2(Xu,Yu,Ustaggered', X, Y, 'cubic');
%V = interp2(Xv,Yv,Vstaggered', X, Y, 'cubic');
%P = interp2(Xp,Yp,Pstaggered', X, Y, 'cubic');


U = sol(uInd)';
V = sol(vInd)';
P = sol(pInd)';

if(toGraph)
        
    figure(1)
    surf(X,Y,P);
    title('pressure');
    xlabel('x');
    ylabel('y');

    figure(2)
    surf(X,Y,U);
    title('u component');
    xlabel('x');
    ylabel('y');

    figure(3)
    surf(X,Y,V);
    title('v component');
    xlabel('x');
    ylabel('y');

    spacing = 4;
    XTemp = X(:,1:spacing:end);
    YTemp = Y(:,1:spacing:end);
    UTemp = U(:,1:spacing:end);
    VTemp = V(:,1:spacing:end);
    
    figure(4);
    quiver(XTemp,YTemp,UTemp,VTemp);
    streamline(X,Y,U,V,0.1,0.5);
end

%return rows / columns that contain all interpolated points
%numYCells = numYCells - 1;
%numXCells = numXCells - 1;

end


