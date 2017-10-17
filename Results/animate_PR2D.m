function [ ] = animate_PR2D( x_traj, u_traj )
    
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
    % Get actual values from symbolics if necessary
    x_traj(1,:) = full(x_traj(1,:));
    x_traj(2,:) = full(x_traj(2,:));
    x_traj(3,:) = full(x_traj(3,:));
    x_traj(4,:) = full(x_traj(4,:));
    x_traj(5,:) = full(x_traj(5,:));
    
    % Rescale
    kilo = 10^3;
    micro = 10^-6;
    x_traj(1,:) = R + kilo * x_traj(1,:);
    x_traj(2,:) = pi/2 + micro * (x_traj(2,:)-pi/2); % subtract pi/2 b.c. initial state
    x_traj(3,:) = kilo * x_traj(3,:);
    x_traj(4,:) = micro * x_traj(4,:);
    
    % Extract trajectory components
    N = length(x_traj(1,:));
    [x, u] = PR2D_pol2cart(x_traj, u_traj);
    pos = x(1:2,:);
    vel = [x(3:4,:); zeros(1,N)];
    mass = x(5,:);
    
    %% Prepare plotting
    % Compute points for target orbit
    angle = linspace(0, 2*pi, 1000);
    r_tar = R + kilo * h_T;
    X_tar = r_tar * cos(angle);
    Y_tar = r_tar * sin(angle);

    % Compute animation timescale
    timeScale = 7 / T;

    % Get number of steps and step width
    N = length(x_traj(1,:));
    DT = T/N;
    
    % Figure properties
    f1 = figure(10);
    clf
    % Set axis background (for demonstration purposes)
    whitebg(f1,[0 0 0]);
    
    %% Animation
    for i=1:N-10
        % Current state
        pos_k = pos(:,i);       % Current position
        pos_mag = norm(pos_k);        % Distance from Body
        vel_k = vel(:,i);       % Current velocity
        
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
        hold on
        % Moon
        rectangle('Position', [-R -R 2*R 2*R], ...
                  'Curvature', [1 1], ...
                  'FaceColor', [0.7 0.7 0.7], ...
                  'LineWidth', 1, ...
                  'EdgeColor', [0.8 0.8 0.8]);
        % Current position
        plot(pos_k(1), pos_k(2), 'x', 'LineWidth', 5);
        % Past trajectory
        plot(pos_traj(1,:), pos_traj(2,:), 'w.', 'LineWidth', 10);
        % Osculating orbit
        plot(X, Y, 'g--');                             
        % Target orbit
        plot(X_tar, Y_tar, 'r');
        % Plot force arrow
        scale = 10^5;
        if i < N
            u_k = u(:,i);
            xArrow = [pos_k(1); pos_k(1) + scale * u_k(1)];
            yArrow = [pos_k(2); pos_k(2) + scale * u_k(2)];
            line(xArrow, yArrow);
        end
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

