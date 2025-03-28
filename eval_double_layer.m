function u = eval_double_layer(eval_pts, Gamma, phi)
%EVAL_DOUBLE_LAYER Evaluates the double-layer potential at given points.
%
% This function computes the value of the double-layer potential at specified
% evaluation points eval_pts using the given boundary parameterization and density phi.
%
% INPUTS:
%   eval_pts - A 2 x m array of evaluation points where the potential is computed.
%              The first row contains x-coordinates, the second row contains y-coordinates.
%   Gamma    - A 2x2 cell array containing function handles for the boundary:
%              Gamma{1,1} : u(t)  -> x-coordinates of boundary points
%              Gamma{1,2} : v(t)  -> y-coordinates of boundary points
%              Gamma{2,1} : u'(t) -> x-component of tangent vector
%              Gamma{2,2} : v'(t) -> y-component of tangent vector
%   phi      - A column vector of length n representing the boundary density function.
%
% OUTPUT:
%   u        - A column vector of size m containing the evaluated double-layer potential.

    n = length(phi);  % Number of boundary discretization points
    
    % Extract function handles for boundary parameterization
    u = Gamma{1,1}; v = Gamma{1,2};
    u_prime = Gamma{2,1}; v_prime = Gamma{2,2};
    
    t = linspace(0, 2 * pi, n+1);  % Parameter t mapped to [0, 2pi]
    t = t(1:end-1);  % Remove duplicate endpoint for periodicity
    t = t';  % Column vector
    dt = t(2) - t(1);  % Uniform spacing

    % Compute boundary points and their derivatives
    x_bdry = u(t);  % x-coordinates of boundary points
    y_bdry = v(t);  % y-coordinates of boundary points
    T_x = u_prime(t);  % x-component of tangent vector
    T_y = v_prime(t);  % y-component of tangent vector
    
    % Compute tangent vector norm |r'(t)|
    T_norm = sqrt(T_x.^2 + T_y.^2);

    % Compute outward normal components
    nu_x = -T_y ./ T_norm;  % Outward normal x-component
    nu_y = T_x ./ T_norm;   % Outward normal y-component

    % Initialize the output array for evaluated values
    u = zeros(size(eval_pts, 2), 1);  % One value for each evaluation point
    
    % Loop over all evaluation points
    for i = 1:size(eval_pts, 2)
        % Extract the current evaluation point
        x_pt = eval_pts(1, i);  % x-coordinate of evaluation point
        y_pt = eval_pts(2, i);  % y-coordinate of evaluation point

        % Compute the integral sum for this evaluation point
        u_val = 0;
        for j = 1:n
            % Compute displacement vector X - Y
            dx = x_pt - x_bdry(j);
            dy = y_pt - y_bdry(j);

            % Compute squared distance |X - Y|^2
            dist_sq = dx^2 + dy^2;

            % Compute normal derivative contribution
            normal_deriv = (dx * nu_x(j) + dy * nu_y(j)) / dist_sq;

            % Compute weighted sum for integral approximation
            u_val = u_val + normal_deriv * phi(j) * T_norm(j) * dt;
        end

        % Store the computed value of u at this point
        u(i) = u_val / (2 * pi);  % Apply normalization factor
    end
end
