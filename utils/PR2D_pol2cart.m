function [ xs_cart, us_cart ] = PR2D_pol2cart( xs_pol, us_pol )
% Takes the trajectory of the PointRocket2D model and converts it
% from polar coordinates to cartesian coordinates
    %% Load model parameters
    run ODEs/PointRocket2D.m
    
    %% Convert states
    % Extract trajectory components and scale to m, rad, kg, s
    rs        = xs_pol(1,:);% * scale(1);
    thetas    = xs_pol(2,:);% * scale(2);
    rDots     = xs_pol(3,:);% * scale(3);
    thetaDots = xs_pol(4,:);% * scale(4);
    ms        = xs_pol(5,:);% * scale(5);
    
    urs       = us_pol(1,:);
    uthetas   = us_pol(2,:);
    
    % Get number of samples
    N = length(rs);
    
    % -- Get position in cartesian frame --
    [posX, posY] = pol2cart(thetas, rs);
    ps = [posX; posY; zeros(1,N)];
    
    % -- Get velocity in cartesian frame --
    % v = v_rel + thetaDots x poss
    vs_rel = [rDots; zeros(2,N)];
    
    % Rewrite angular velocity as rotation around z-axis
    thetaDots = [zeros(2,N); thetaDots]; 
    
    % Compute velocity component of rotation
    vs_rot = cross(thetaDots, ps);
    
    % Compute velocities
    vs = vs_rel + vs_rot;
    
    % -- Return state trajectory in cartesian coordinates
    xs_cart = [ps(1:2,:); vs(1:2,:); ms];
    
    %% Convert controls
    % Convert to cartesian coordinates (note the corotating reference frame)
    us_cart = [urs; uthetas];
    %us_cart = [uthetas; urs];

    % Rotate into inertial frame
    for i = 1:size(us_pol,2)
        theta = thetas(i);
        Rotation = [cos(theta), -sin(theta); ...
                    sin(theta),  cos(theta)];
        
        us_cart(:,i) = Rotation * us_cart(:,i);
    end
end

