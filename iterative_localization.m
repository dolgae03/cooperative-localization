function [h_final, r_final] = iterative_localization(h, r, epsilon, gamma_initial, D, delta, Q1, Q2)
    % h: Initial estimates of the positions.
    % r: Initial estimates of the range.
    % rij: Matrix of range measurements between points.
    % C: Covariance matrix.
    % b: Vector b from the optimization problem.
    % gamma_initial: Initial gamma value.
    % D: Prior error range of the GNSS information.
    % delta: Convergence tolerance.
    % max_iterations: Maximum number of iterations.


    %% Step 1: Initalize value
    % Initialize variables
    k = 0;
    u_hat_k_1 = [h', r']'; % Initial u_k = [h', r']'
    
    b = [h', r']';
    n = height(h) / 2;
    m = height(r);
    C = [inv(Q1) zeros(2*n,m);
         zeros(m,2*n) inv(Q2)];
    
    
    while true
        %% Step 2: Update Ak-1 and Pk-1
        gamma_k = gamma_initial;
        k = k + 1;
        A_k_1 = update_A(n, m, u_hat_k_1, epsilon);
        P_k_1 = eye(size(A_k_1, 2)) - A_k_1' * inv(A_k_1 * A_k_1') * A_k_1;
        
        %% Step 3: Solve the approximate problem
        u_bar_k_1 = pinv(P_k_1 * C * P_k_1) * C * b;


        for threshold = 1:15
            %% Step 4: Calculate the update of the kth iteration
            u_hat_k = u_hat_k_1 + gamma_k * (u_bar_k_1 - u_hat_k_1);
            
            %% Step 5: Inspection step I
            x_estimated = u_hat_k(1:length(h))';
            
            diff_x = reshape(x_estimated' - h, [], 2);
            % max(vecnorm(diff_x, 2, 2))
            if max(vecnorm(diff_x, 2, 2)) > D
                gamma_k = gamma_k / 2;
                continue; % go to Step 4
            end
            
            %% Step 6: Correction step
            new_d = correct_d(x_estimated, epsilon)';
            u_hat_k(length(h)+1:end) = new_d;
            
            %% Step 7: Inspection step II
            objective_k = u_hat_k' * C * u_hat_k - 2 * b' * C * u_hat_k;
            objective_k_1 = u_hat_k_1' * C * u_hat_k_1 - 2 * b' * C * u_hat_k_1;
            if objective_k > objective_k_1
                gamma_k = gamma_k / 2;
                continue; % go to Step 4
            end

            break
        end
  
        %% Step 8: Convergence check
        relative_error = norm(u_hat_k - u_hat_k_1) / norm(u_hat_k);
        if relative_error <= delta
            break;
        end

        u_hat_k_1 = u_hat_k;
    end
    
    u_final = u_hat_k; % Final output

    h_final = u_hat_k(1:length(h));
    r_final = u_hat_k(length(h)+1:end);
end

function C = get_C_matrix(n, m, i, j, epsilon)
    % epsilon is size x 2 matrix describe edge between,.
    % Update A_k_minus_1 based on the definitions provided.
    C = zeros(2*n + m, 2*n + m);

    C(2*i-1:2*i, 2*i-1:2*i) = eye(2);
    C(2*j-1:2*j, 2*j-1:2*j) = eye(2);
    C(2*i-1:2*i, 2*j-1:2*j) = -eye(2);
    C(2*j-1:2*j, 2*i-1:2*i) = -eye(2);
   
    sum_nl = sum(epsilon(:, 1) < i);

    C(2*n + sum_nl + j, 2*n + sum_nl + j) = -1;
end


function A_k = update_A(n, m, u_k, epsilon)
    A_k = zeros(height(epsilon), height(u_k));
    

    for idx = 1:height(epsilon)
        i = epsilon(idx, 1);
        j = epsilon(idx, 2);
        
        Cij = get_C_matrix(n, m, i, j, epsilon);

        A_k(idx, :) = u_k' * Cij;
    end
end

function d = correct_d(x_estimated, epsilon)
    d = zeros(height(epsilon), 1);
    for idx=1:height(epsilon)
        i = epsilon(idx, 1);
        j = epsilon(idx, 2);

        d(idx, 1) = norm(x_estimated(2*i-1 : 2*i) - x_estimated(2*j-1 : 2*j));
    end
end