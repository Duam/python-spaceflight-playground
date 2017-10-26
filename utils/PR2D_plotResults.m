function [] = PR2D_plotResults(sol)
% Plot the trajectory
    
    coordSys = sol.param.coordSys;

    % Load parameters
    % Current working directory must be "ControlledRocket"
    if strcmpi(coordSys, 'pol')
        run ODEs/PointRocket2D.m
    elseif strcmpi(coordSys, 'cart')
        run ODEs/PointRocket2D_cart.m
    end

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
    angVel_T = sol.param.thetaDotT;

    % Time axis for plotting
    tAxis = 0:DT:T-DT;
    
    if strcmpi(coordSys, 'pol')
        % Radius of planet
        r_0 =      zeros(1,length(tAxis));
        r_T = hT * ones(1,length(tAxis));
        % Angle
        theta_0 = ones(1,length(tAxis));
        % Angular velocity limits
        thetaDot_0 = zeros(1,length(tAxis));
        thetaDot_T = angVel_T * ones(1,length(tAxis));
        % Mass limits
        m_0 = m0 * ones(1,length(tAxis));
        m_T = m_e * ones(1,length(tAxis));

        figure(1);
        clf
        % Plot radius
        subplot(3,2,1);
        hold on
        plot(tAxis, r_0, '--r');
        plot(tAxis, r_T, '--r');
        plot(tAxis, X(1,:));
        ylabel('$r$', 'interpreter', 'latex');
        grid on
        % Plot angle
        subplot(3,2,2);
        hold on
        plot(tAxis, theta_0, '--r');
        plot(tAxis, X(2,:));
        ylabel('$\theta$', 'interpreter', 'latex');
        grid on
        % Plot radial velocity
        subplot(3,2,3);
        hold on
        plot(tAxis, zeros(1,N), '--r');
        plot(tAxis, X(3,:));
        ylabel('$\dot{r}$', 'interpreter', 'latex');
        grid on
        % Plot angular velocity
        subplot(3,2,4);
        hold on
        plot(tAxis, thetaDot_0, '--r');
        plot(tAxis, thetaDot_T, '--r');
        plot(tAxis, X(4,:));
        ylabel('$\dot{\theta}$', 'interpreter', 'latex');
        grid on
        % Plot mass
        subplot(3,2,[5 6]);
        hold on
        plot(tAxis, m_0, '--r');
        %plot(tAxis, m_T, '--r');
        plot(tAxis, X(5,:));
        ylabel('$m$', 'interpreter', 'latex');
        grid on

        %% Plot controls
        figure(2);
        clf
        % Plot radius
        subplot(2,1,1);
        hold on
        stairs(tAxis(1:end), U(1,:));
        ylabel('$u_r$', 'interpreter', 'latex');
        grid on
        % Plot angle
        subplot(2,1,2);
        hold on
        stairs(tAxis(1:end), U(2,:));
        ylabel('$u_r$', 'interpreter', 'latex');
        grid on
        
    elseif strcmpi(coordSys, 'cart')
        figure(1);
        clf
        % Plot position
        subplot(2,2,1);
        hold on
        plot(x0(1), x0(2), 'o');
        plot(X(1,:), X(2,:));
        xlabel('$p_x$', 'interpreter', 'latex');
        ylabel('$p_y$', 'interpreter', 'latex');
        grid on
        % Plot velocity
        subplot(2,2,2);
        hold on
        plot(x0(3), x0(4), 'o');
        plot(X(3,:), X(4,:));
        xlabel('$v_x$', 'interpreter', 'latex');
        ylabel('$v_y$', 'interpreter', 'latex');
        grid on
        % Plot mass
        subplot(2,2,[3 4]);
        hold on
        plot(0, x0(5), 'o');
        plot(tAxis, X(5,:));
        ylabel('$m$', 'interpreter', 'latex');
        grid on

        %% Plot controls
        figure(2);
        clf
        % Plot radius
        subplot(2,1,1);
        hold on
        stairs(tAxis(1:end), U(1,:));
        ylabel('$u_x$', 'interpreter', 'latex');
        grid on
        % Plot angle
        subplot(2,1,2);
        hold on
        stairs(tAxis(1:end), U(2,:));
        ylabel('$u_y$', 'interpreter', 'latex');
        grid on
    end
end