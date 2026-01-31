%% 1. Initialization and Parameters
clear variables; clc;

% System Physical Parameters
m1 = 1.2;       % Mass of cart
m2 = 0.6;       % Mass of pendulum
l  = 1.5;       % Length of pendulum
g  = 9.8067;    % Gravity

% Controller / Planner Settings
ctrl_intervals = 100;       % Number of control intervals (resolution)
ctrl_const     = [-2, 2];   % Input constraints [min, max]
createPath     = false;     % Flag: true = generate new path, false = load existing

%% 2. Path Generation or Loading
if createPath
    % --- A. Setup Environment ---
    % Define paths to dependencies
    yop_path    = '/home/filipe/Documents/Research/smpc-with-risk-allocation/Example 2/model_in_yop/yop';
    casadi_path = '/home/filipe/Documents/casadi-3.7.2';
    
    addpath(genpath(yop_path));
    addpath(genpath(casadi_path));

    % --- B. Generate Optimal Path ---
    % Calls the external Yop/CasADi function
    path = pathplanner(m1, m2, l, g, ctrl_intervals, ctrl_const);
    
    % Clean up path after generation to avoid conflicts later
    restoredefaultpath; 

    % --- C. Data Extraction & Resampling ---
    % Construct matrix from planner output
    % Row 1: Independent variable (Time)
    % Row 2-4: States (Velocity, Angle, AngVelocity) 
    % Note: State(1) (Position) is used separately as the interpolation base
    x_path_original = [path.Independent;
                       path.State(2,:);
                       path.State(3,:);
                       path.State(4,:)];
    
    % Use Position (State 1) as the reference for resampling
    time_irr = path.State(1,:); 
    
    % Calculate sampling frequency based on total distance/intervals
    fs = ctrl_intervals / time_irr(end);
    Ts = time_irr(end) / ctrl_intervals;

    % Resample trajectory to fixed grid
    % 1. Resample Time vector
    [~, time] = resample(x_path_original(1,:), time_irr, fs);
    
    % 2. Resample States (Time, Velocity, Angle, AngVel)
    x_path = [resample(x_path_original(1,:), time_irr, fs);
              resample(x_path_original(2,:), time_irr, fs);
              resample(x_path_original(3,:), time_irr, fs);
              resample(x_path_original(4,:), time_irr, fs)];
          
    % 3. Resample Inputs
    u_path = resample(path.Input, time_irr, fs);

else
    % --- D. Load Pre-calculated Path ---
    % Loads: x_path, u_path, Ts, time, and potentially original parameters
    load path_stored.mat
end