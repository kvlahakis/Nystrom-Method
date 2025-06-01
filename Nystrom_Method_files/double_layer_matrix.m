function K = double_layer_matrix(Gamma, n)

    % DOUBLE_LAYER_MATRIX Constructs the double-layer potential matrix.
    % 
    % This function computes the matrix K corresponding to the double-layer 
    % potential integral operator for a given boundary parameterization. 
    % The double-layer potential is used to solve boundary integral equations 
    % for Laplace's equation in a 2D domain.
    %
    % INPUTS:
    %   Gamma - A 3x2 cell array containing function handles for the boundary:
    %           Gamma{1,1} : u(t)  -> x-coordinates of boundary points
    %           Gamma{1,2} : v(t)  -> y-coordinates of boundary points
    %           Gamma{2,1} : u'(t) -> x-component of tangent vector
    %           Gamma{2,2} : v'(t) -> y-component of tangent vector
    %           Gamma{3,1} : u''(t) -> x-component of second derivative
    %           Gamma{3,2} : v''(t) -> y-component of second derivative
    %
    %   n     - Number of discretization points along the boundary
    %
    % OUTPUT:
    %   K     - n x n matrix representing the discretized double-layer operator
    
    % Extract function handles for u(t), v(t), u'(t), v'(t)
    u = Gamma{1,1}; v = Gamma{1,2};
    u_prime = Gamma{2,1}; v_prime = Gamma{2,2};
    u_dprime = Gamma{3,1}; v_dprime = Gamma{3,2};
    
    t = linspace(0, 2 * pi, n+1);  % Parameter t mapped to [0, 2pi]
    t = t(1:end-1);
    t = t';    
    dt = t(2) - t(1);  % Uniform spacing

    % Compute boundary points r(t) and derivatives r'(t)
    x = u(t);  % x-coordinates of boundary points
    y = v(t);  % y-coordinates of boundary points
    T_x = u_prime(t);  % x-component of tangent vector
    T_y = v_prime(t);  % y-component of tangent vector
    
    u_ddot = u_dprime(t); % x-component of second derivative of boundary point
    v_ddot = v_dprime(t); % y-component of second derivative of boundary point
    
    % Compute normalizing factor |r'(t)|
    T_norm = sqrt(T_x.^2 + T_y.^2);

    % Compute outward normal components
    nu_x = -T_y ./ T_norm;  % Outward normal x-component
    nu_y = T_x ./ T_norm;   % Outward normal y-component

    % Initialize the matrix K
    K = zeros(n, n);

    % Compute the matrix elements
    for i = 1:n % fix a target point X = (u(s), v(s))
        for j = 1:n % fix an integration node Y = (u(t), v(t))
            if i ~= j
                % Compute X - Y
                dx = x(i) - x(j);
                dy = y(i) - y(j);
                
                % Compute distance squared |X - Y|^2
                dist_sq = dx^2 + dy^2;

                % Compute (X - Y) . nu_y (dot product)
                dot_prod = dx * nu_x(j) + dy * nu_y(j);

                % Compute kernel entry
                K(i, j) = dot_prod / (dist_sq * 2 * pi) * T_norm(j) * dt;
            end
        end
    end

    for i = 1:n
        K(i, i) = -1/(4*pi * T_norm(i)^2) * (T_y(i) * u_ddot(i) - T_x(i) * v_ddot(i)) * dt;
    end
end