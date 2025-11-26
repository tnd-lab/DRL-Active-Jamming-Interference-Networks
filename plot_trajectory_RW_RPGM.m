close all
clear all

%% Jammer
pos_jammer = load('data_pos_RWP.csv');

vX = pos_jammer(:,1); 
vY = pos_jammer(:,2); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
jammer = plot(vX, vY, 'rd:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;

%%
pos_nodes = load('data_pos_RPGM.csv');

%% Controller Tx0
vX = pos_nodes(:,1); 
vY = pos_nodes(:,2); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Tx0 = plot(vX, vY, 'gs:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;

%% Tx1
vX = pos_nodes(:,5); 
vY = pos_nodes(:,6); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Tx1 = plot(vX, vY, 'b^:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;

%% Rx1
vX = pos_nodes(:,7); 
vY = pos_nodes(:,8); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Rx1 = plot(vX, vY, 'bo:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;

%% Tx2
vX = pos_nodes(:,9); 
vY = pos_nodes(:,10); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Tx2 = plot(vX, vY, 'b^:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;
grid on

%% Rx2
vX = pos_nodes(:,11); 
vY = pos_nodes(:,12); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Rx2 = plot(vX, vY, 'bo:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;

%% Tx3
vX = pos_nodes(:,13); 
vY = pos_nodes(:,14); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Tx3 = plot(vX, vY, 'b^:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;
grid on

%% Rx3

vX = pos_nodes(:,15); 
vY = pos_nodes(:,16); 
rMag = 0.5;
% Length of vector
lenTime = length(vX);
% Indices of tails of arrows
vSelect0 = 1:(lenTime-1);
% Indices of tails of arrows
vSelect1 = vSelect0 + 1;
% X coordinates of tails of arrows
vXQ0 = vX(vSelect0, 1);
% Y coordinates of tails of arrows
vYQ0 = vY(vSelect0, 1);
% X coordinates of heads of arrows
vXQ1 = vX(vSelect1, 1);
% Y coordinates of heads of arrows
vYQ1 = vY(vSelect1, 1);
% vector difference between heads & tails
vPx = (vXQ1 - vXQ0) * rMag;
vPy = (vYQ1 - vYQ0) * rMag;
% make plot 
Rx3 = plot(vX, vY, 'bo:', 'linewidth', 0.5, 'markersize', 5); hold on;
% add arrows 
quiver(vXQ0,vYQ0, vPx, vPy, 0, 'r'); hold on;


%%
grid on
xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')
legend([jammer, Tx0, Tx1, Rx1], ...
    {'Jammer', 'Hub', ...
    '$\mathrm{T}_n$', '$\mathrm{R}_n$'},...
    'location', 'best',...
    'interpreter', 'latex')
set(gca,'FontSize',17)
set(gca, 'LooseInset', get(gca, 'TightInset')) %remove plot padding
axis([0 500 0 500])
xticks = get(gca,'XTickLabel');
set(gca,'XTickLabel',xticks,'FontName','Times','fontsize',15)