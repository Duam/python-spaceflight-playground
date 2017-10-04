function [ ] = animate_PR2D( x_traj )
    
    %% Load ODE and parameters
    % Current working directory must be "ControlledRocket"
    run ODEs/PointRocket2D.m
    
    % Simulation parameters
    T = 600;
    N = 100;
    DT = T/N;

    % Target altitude
    h_T = 20;
        
    %% Get trajectory in euclidian coordinate frame (remember scaling)
    kilo = 10^3;
    micro = 10^-6;
    x_traj(1,:) = kilo * (R + x_traj(1,:));
    x_traj(2,:) = pi/2 + micro * (x_traj(2,:)-pi/2); % subtract pi/2 b.c. initial state
    x_traj(3,:) = kilo * x_traj(3,:);
    x_traj(4,:) = micro * x_traj(4,:);
    
    % Extract trajectory components (account for scaling)
    N = length(x_traj(1,:));
    [posX, posY] = pol2cart(x_traj(2,:),x_traj(1,:));
    pos = [posX; posY];
    
    % Compute linear velocity
    % v = v_rel + w x r
    v_rel = [x_traj(3,:); zeros(2,N)];
    
    % Rewrite angular velocity as rotation around z-axis
    omega = [zeros(2,N); x_traj(4,:)];
    
    % Compute velocity component of rotation
    v_rot = cross(omega, [pos;zeros(1,N)]);
    
    % Compute velocity
    vel = v_rel + v_rot;
   
    % Get mass
    mass = x_traj(5,:);

    %% Prepare plotting
    % Compute points for target orbit
    angle = linspace(0, 2*pi, 1000);
    r_tar = kilo * (R + h_T);
    X_tar = r_tar * cos(angle);
    Y_tar = r_tar * sin(angle);

    % Compute animation timescale
    timeScale = 1 / T;

    % Get number of steps and step width
    N = length(x_traj(1,:));
    DT = T/N;
    
    % Figure properties
    f1 = figure(10);
    clf
    % Set axis background (for demonstration purposes)
    whitebg(f1,[0 0 0]);
    
    %% Animation
    for i=1:N
        % Current state
        pos_k = full(pos(:,i));       % Current position
        pos_mag = norm(pos_k);        % Distance from Body
        vel_k = full(vel(:,i));       % Current velocity
        
        % All states up until now
        pos_traj = [pos(:,1:i), inf(2,N-i)];
        
        % Compute Kepler elements
        h = cross([pos_k; 0],vel_k);
        mu = G * (M + mass(i));
        e = cross(vel_k,h)/mu - [pos_k; 0]/pos_mag;
        p = norm(h)^2 / mu;     % Something rectum
        a = p / (1-norm(e)^2);  % Length of big semi-axis
        b = sqrt(p*a);          % Length of small semi-axis
        
        % Excentricity points towards periapsis. Used for rotating ellipse
        e_thet = atan2(e(2), e(1)); % Angle of excentricity

        % Elements of osculating orbit
        pos_pe  = e/norm(e) * p/(1+norm(e));  % Periapsis
        pos_ap  = -e/norm(e) * p/(1-norm(e)); % Apoapsis
        pos_mid = (pos_ap + pos_pe) / 2;      % Midpoint of ellipse
        pos_f1  = [0;0;0];                    % Focus point 1 (body)
        pos_f2  = pos_f1 + pos_mid * 2;       % Focus point 2
        
        % Compute osculating orbit
        beta = e_thet;
        x = pos_mid(1);
        y = pos_mid(2);
        [X,Y] = computeEllipse(x, y, a, b, beta);
       
        %% Plotting
        % Make figure or use existing one
        figure(f1);
        clf;
        
        % Draw elements
        hold on
        rectangle('Position', kilo*[-R -R 2*R 2*R], ...
                  'Curvature', [1 1], ...
                  'FaceColor', [0.7 0.7 0.7], ...
                  'LineWidth', 1, ...
                  'EdgeColor', [0.8 0.8 0.8]);
        plot(pos_k(1), pos_k(2), 'x', 'LineWidth', 5); % Current position
        plot(pos_traj(1,:), pos_traj(2,:), 'w.', 'LineWidth', 10);       % Past trajectory
        plot(X, Y, 'g--');                             % Osculating orbit
        % Plot target orbit
        plot(X_tar, Y_tar, 'r');                       % Target Orbit
        hold off
        
        % Figure properties
        % Follow rocket
        xMin = pos_k(1) - 1 * 10^4; %, x_traj(1,:)]);
        xMax = pos_k(1) + 1 * 10^4; %, x_traj(1,:)]);
        yMin = pos_k(2) - 1 * 10^4; %, x_traj(2,:)]);
        yMax = pos_k(2) + 1 * 10^4; %, x_traj(2,:)]);
%         Fixed view
%         xMin = -r_tar/2 - 5 * 10^4;
%         xMax = 0    + 5 * 10^4;
%         yMin = r_tar-20000    - 5 * 10^4;
%         yMax = r_tar   + 5 * 10^4;
%         
        axis([xMin xMax yMin yMax]);
%         axis equal
        title('Trajectory');
        xlabel('x');
        ylabel('y');
        grid on
        drawnow
        pause(DT*timeScale);
        
        
    end
    
end

