clear all
close all
clc

%2D Stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the L2 error with the actual analytical solution for g = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num = 40;
testPoints = num-20:1:num;

L2EP = zeros(length(testPoints),1);
L2EU = zeros(length(testPoints),1);
L2EV = zeros(length(testPoints),1);
d = zeros(length(testPoints),1);

width = 20;
height = 20;
L = 5;
R = 5;
g = 0;
mu = 2.0;
p0 = 200;


% Compute the L2 error
%2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
for i = 1:length(testPoints)
    
    %function [ P U V X Y numYCells numXCells d] = StokesStaggered(g, numYCells, p0, mu, toGraph, height, width )
    [ p u v xu yu xv yv xp yp numYCells numXCells delta] = StokesStaggered(g, testPoints(i), p0, mu, g, height, width, R, L);
    [P U V] = pTest(xu,yu,xv,yv,xp,yp,R,L);
    
    d(i) = delta;
    
    L2EP(i) = sum(sum((p - P).^2));
    L2EU(i) = sum(sum((u - U).^2));
    L2EV(i) = sum(sum((v - V).^2));

    L2EP(i) = sqrt(L2EP(i) / prod(size(xp)));
    L2EU(i) = sqrt(L2EU(i) / prod(size(xu)));
    L2EV(i) = sqrt(L2EV(i) / prod(size(xv)));
end

d = d(4:end);
L2EP = L2EP(4:end);
L2EU = L2EU(4:end);
L2EV = L2EV(4:end);

figure()
subplot(2,2,1);
loglog(d,L2EP,'-',d, d.^2 *  L2EP(end) / d(end)^2,'--');
%loglog(dx,L2EP,'-',dx,dx.^2,'--');
title('Error for P (Pressure)');

subplot(2,2,2);
loglog(d,L2EU,'-',d,d.^2 *  L2EU(end) / d(end)^2,'--');
title('Error for U (Horizontal Velocity)');

subplot(2,2,3);
loglog(d,L2EV,'-',d,d.^2 *  L2EV(end) / d(end)^2,'--');
title('Error for V (Vertical Velocity)');

figure()
surf(xu,yu,u);
title('u');
xlabel('x');
ylabel('y');

figure()
surf(xv,yv,v);
title('v');
xlabel('x');
ylabel('y');

figure()
surf(xp,yp,p);
title('p');
xlabel('x');
ylabel('y');

%%%%%%%%%%%%%benchmark%%%%%%%%%%%%%%%%
REPS = 10;
num = [25 50 75 100 200];
tic;
for i = 1:length(num)
    num(i)
    for j = 1:REPS
        [ p u v x y] = StokesStaggered(g, num(i), p0, mu, g, height, width, R, L);
    end
    averageTime(i) = toc / REPS;
end

averageTime

figure()
plot(num,averageTime);
title('Execution Time vs. Number of Discretization Points');
xlabel('M');
ylabel('Execution Time (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





