clc; close all;
%%Init
syms x1 x2;
w = 0.9;    %inertial coefficient
c1 = 1.9;   %cognitive coefficient
c2 = 1.9;   %social coeffiecient
numPos = 100;   %number of positions
epochs = 50;    %number of iterations
%define objective function
f = @(x1, x2) ((x2-x1).^4 + 12.*x1.*x2 - x1 + x2 - 3);
%empty cell arrays to store x_positions, p_positions, velocities and
%function evaluations and gbest, worst, average values
X_pos = cell(epochs,1); 
p_pos = cell(epochs,1);
vel = cell(epochs,1);
fvalx = cell(epochs,1);
gBest = cell(epochs,1);
bestVal = zeros(epochs,1);
worsVal = zeros(epochs,1);
meanVal = zeros(epochs,1);
X_pos{1} = rand(numPos,2)*2.0-1;    %Initialize random positions of points
vel{1} = rand(numPos,2);    %initialize random velocities
%initialize random r and s
r = rand(numPos,2);
s = rand(numPos,2);
fvalx{1} = f(X_pos{1}(:,1),X_pos{1}(:,2));  %evaluate objfunc values
p_pos{1} = X_pos{1};    %set p_pos 
[~,ind] = min(fvalx{1});    %get minimum func evaluation to set gbest
gBest{1} = X_pos{1}(ind,:); 
bestVal(1) = f(gBest{1}(:,1), gBest{1}(:,2));
worsVal(1) = max(fvalx{1});
meanVal(1) = mean(fvalx{1});
%% Run for k iterations
for k = 1:epochs-1
    vel{k+1} = (w.*vel{k})+(c1.*r.*(p_pos{k}-X_pos{k}))...
        +(c2.*s.*(gBest{k}-X_pos{k}));  %update velocities
    X_pos{k+1} = X_pos{k}+vel{k+1}; %update positions 
    fvalx{k+1} = f(X_pos{k+1}(:,1),X_pos{k+1}(:,2));    %evaluate function
    p_pos{k+1} = p_pos{k};  %update pbest
    for i = 1:numPos
        if fvalx{k+1}(i)<f(p_pos{k}(i,1),p_pos{k}(i,2))
            p_pos{k+1}(i,:) = X_pos{k+1}(i,:);  %update pbest 
        end
    end
    [~,ind] = min(fvalx{k+1});  %min func eval for gbest
    gBest{k+1} = X_pos{k+1}(ind,:); 
    bestVal(k+1) = f(gBest{k+1}(:,1), gBest{1}(:,2));
    worsVal(k+1) = max(fvalx{k+1});
    meanVal(k+1) = mean(fvalx{k+1});
end
%% Plotting
X = linspace(-1,1,1000);
Y = X;
[Xpt, Ypt] = meshgrid(X,Y);
Z = f(Xpt,Ypt);
figure;
contour(Xpt,Ypt, Z, 30); hold on;
plot(gBest{end}(1),gBest{end}(2), 'o');
figure;
xax = 1:epochs; %x axis
% worst, best and average value of the functions at every point
plot(xax, worsVal, 'r-', xax, meanVal, 'b-', xax, bestVal, 'g-');
legend('Worst Value', 'Mean Value', 'Best Value');

%% Scatter plot every iteration to visualize the convergence
% for i = 1:length(X_pos)
% figure;
% scatter(X_pos{i}(:,1), X_pos{i}(:,2));
% xlim([-1 1]); ylim([-1 1]);
% name = horzcat('it',num2str(i), '.png');
% saveas(gcf, name);
% end