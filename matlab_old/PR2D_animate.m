function [ ] = PR2D_animate( sol )
% Animate the trajectory    

    % Coordinate system information
    coordSys = sol.param.coordSys;

    % Trajectories
    x0 = sol.x0;
    X = sol.X;
    U = sol.U;
    
    % Simulation parameters
    T = sol.param.T;
    N = sol.param.N;
    DT = T/N;

    % Target altitude
    hT = sol.param.hT;
    % Orbital angular velocity in microradians
    angVelT = sol.param.thetaDotT;
    
    % Load parameters & convert
    % Current working directory must be "ControlledRocket"
    if strcmpi(coordSys, 'pol')
        run ODEs/PointRocket2D.m
        [x0, u0] = PR2D_pol2cart(x0, [0;0]);
        [X, U] = PR2D_pol2cart(X, U); % Convert to cartesian
    elseif strcmpi(coordSys, 'cart')
        run ODEs/PointRocket2D_cart.m
    end

    % Split trajectory vector
    ps = X(1:2,:);
    vs = X(3:4,:);
    ms = X(5,:);
    
    % Prepare plotting
    % Compute points for target orbit
    angle = linspace(0, 2*pi, 1000);
    r_orb = R + 10^3 * hT;
    X_orb = r_orb * cos(angle);
    Y_orb = r_orb * sin(angle);

    % Compute animation timescale
    timeScale = 7 / T;
    
    % Figure properties
    f1 = figure(10);
    clf
    % Set axis background (for demonstration purposes)
    whitebg(f1,[0 0 0]);
    
    % Animation
    for i=1:N
        % Current state
        pk = [ps(:,i); 0];       % Current position (3D)
        pk_mag = norm(pk);       % Distance from origin
        vk = [vs(:,i); 0];       % Current velocity (3D)
        mk = ms(i);              % Current mass
        
        % All positions up until now
        ps_past = [ps(:,1:i), inf(2,N-i)];
        
        % Compute Kepler elements
        h = cross(pk,vk);
        mu = G * (M + mk);
        e = cross(vk,h)/mu - pk/pk_mag;
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
%         pos_f2  = pos_f1 + pos_mid * 2;       % Focus point 2
        
        % Compute osculating orbit
        beta = e_thet;
        X_osc = pos_mid(1);
        Y_osc = pos_mid(2);
        [X,Y] = computeEllipse(X_osc, Y_osc, a, b, beta);
       
        
        % Plotting
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
        plot(pk(1), pk(2), 'x', 'LineWidth', 5);
        % Past trajectory
        plot(ps_past(1,:), ps_past(2,:), 'w.', 'LineWidth', 10);
        % Osculating orbit
        plot(X, Y, 'g--');                             
        % Target orbit
        plot(X_orb, Y_orb, 'r');
        % Plot force arrow
        scale = 10^4;
        if i < N
            u_k = U(:,i);
            xArrow = [pk(1); pk(1) + scale * u_k(1)];
            yArrow = [pk(2); pk(2) + scale * u_k(2)];
            line(xArrow, yArrow);
        end
        hold off
        
        % Figure properties
        % Follow rocket
        xMin = pk(1) - 1 * 10^4; %, x_traj(1,:)]);
        xMax = pk(1) + 1 * 10^4; %, x_traj(1,:)]);
        yMin = pk(2) - 1 * 10^4; %, x_traj(2,:)]);
        yMax = pk(2) + 1 * 10^4; %, x_traj(2,:)]);
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

