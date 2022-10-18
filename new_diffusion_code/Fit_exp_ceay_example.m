% Uses fitnlm() to fit a non-linear model (an exponential decay curve, Y = a * exp(-b*x)) through noisy data.
% Requires the Statistics and Machine Learning Toolbox, which is where fitnlm() is contained.
% Initialization steps.
clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 20;
% Create the X coordinates from 0 to 20 every 0.5 units.
X = 0 : 0.5 : 20;
% Define coefficients and function that the X values obey.
a = 10 % Arbitrary sample values I picked.
b = 0.4
Y = a * exp(-X * b); % Get a vector.  No noise in this Y yet.
% Add noise to Y.
Y = Y + 0.5 * randn(1, length(Y));
%--------------------------------------------------------------------------------------------------------------------------------------
% Now we have noisy training data that we can send to fitnlm().
% Plot the noisy initial data.
plot(X, Y, 'b*', 'LineWidth', 2, 'MarkerSize', 15);
grid on;
% Convert X and Y into a table, which is the form fitnlm() likes the input data to be in.
% Note: it doesn't matter if X and Y are row vectors or column vectors since we use (:) to get them into column vectors for the table.
tbl = table(X(:), Y(:));
% Define the model as Y = a * exp(-b*x)
% Note how this "x" of modelfun is related to big X and big Y.
% x((:, 1) is actually X and x(:, 2) is actually Y - the first and second columns of the table.
modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1));  
% Guess values to start with.  Just make your best guess.
aGuessed = 11 % Arbitrary sample values I picked.
bGuessed = 0.5
beta0 = [aGuessed, bGuessed]; % Guess values to start with.  Just make your best guess.
% Now the next line is where the actual model computation is done.
mdl = fitnlm(tbl, modelfun, beta0);
% Now the model creation is done and the coefficients have been determined.
% YAY!!!!
% Extract the coefficient values from the the model object.
% The actual coefficients are in the "Estimate" column of the "Coefficients" table that's part of the mode.
coefficients = mdl.Coefficients{:, 'Estimate'}
% Create smoothed/regressed data using the model:
yFitted = coefficients(1) * exp(-coefficients(2)*X);
% Now we're done and we can plot the smooth model as a red line going through the noisy blue markers.
hold on;
plot(X, yFitted, 'r-', 'LineWidth', 2);
grid on;
title('Exponential Regression with fitnlm()', 'FontSize', fontSize);
xlabel('X', 'FontSize', fontSize);
ylabel('Y', 'FontSize', fontSize);
legendHandle = legend('Noisy Y', 'Fitted Y', 'Location', 'north');
legendHandle.FontSize = 30;
% Place formula text roughly in the middle of the plot.
formulaString = sprintf('Y = %.3f * exp(-%.3f * X)', coefficients(1), coefficients(2))
xl = xlim;
yl = ylim;
xt = xl(1) + abs(xl(2)-xl(1)) * 0.325;
yt = yl(1) + abs(yl(2)-yl(1)) * 0.59;
text(xt, yt, formulaString, 'FontSize', 25, 'FontWeight', 'bold');
% Set up figure properties:
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 