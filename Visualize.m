close all;
clear all;
tic

% Problem Definition
n = 2; N = 3; T = 30; D = 0.4; T_s = 0.1;

% Define our cost matrices
Q = diag([.1 .1 0 0]) * 1.0; % cost associated with the state
Qf = eye(2*n,2*n) * 10.0; % cost associated with the final state
R = eye(n, n) * 0.1; % cost associated with the control effort 
collision_scalar = 10.0;

% Dynamics matrices
A_d = [eye(n, n) T_s*eye(n,n); zeros(n,n) eye(n,n)];
B_d = [T_s^2/2*eye(n,n); T_s*eye(n,n)];

% Control constraints
u_max = 2.0;
u_min = -2.0;

% Initial and final positions
p1init = [-1,0, 0, 0]'; p2init = [-1,0.5, 0, 0]'; p3init = [-1,1, 0, 0]';
p1final = [1,0, 0, 0]'; p2final = [1,-0.5, 0, 0]'; p3final = [1,-1, 0, 0]';

% Initial trajectory
cur_p1 = repmat(p1init,1,T)+(p1final-p1init)*[0:1/(T-1):1];
cur_p2 = repmat(p2init,1,T)+(p2final-p2init)*[0:1/(T-1):1]; 
cur_p3 = repmat(p3init,1,T)+(p3final-p3init)*[0:1/(T-1):1]; 

% Warm start trajectory
[cur_p1, cur_p2, cur_p3, ~] = traj_opt_fcn(p1init(1:n, :), p2init(1:n, :), p3init(1:n, :), p1final(1:n, :), p2final(1:n, :), p3final(1:n, :), cur_p1(1:n, :), cur_p2(1:n, :), cur_p3(1:n, :), 10);
cur_p1(n+1:2*n,:) = [(cur_p1(1:n,2:T)-cur_p1(1:n,1:T-1))/T_s [0;0]]; 
cur_p2(n+1:2*n,:) = [(cur_p2(1:n,2:T)-cur_p2(1:n,1:T-1))/T_s [0;0]]; 
cur_p3(n+1:2*n,:) = [(cur_p3(1:n,2:T)-cur_p3(1:n,1:T-1))/T_s [0;0]]; 

% random
%function [p1_ret, p2_ret, p3_ret, costs] = traj_opt_fcn(p1init, p2init, p3init, p1final, p2final, p3final, cur_p1, cur_p2, cur_p3, iterations)

% cur_p1 = randn(2*n, T); cur_p1(:, 1) = p1init; cur_p1(:, end) = p1final;
% cur_p2 = randn(2*n, T); cur_p2(:, 1) = p2init; cur_p2(:, end) = p2final;
% cur_p3 = randn(2*n, T); cur_p3(:, 1) = p3init; cur_p3(:, end) = p3final;

iterations = 2;

costs = zeros(1, iterations+1);
costs(1) = norm(cur_p1(:, 2:end) - cur_p1(:, 1:end-1), 'fro') + norm(cur_p2(:, 2:end) - cur_p2(:, 1:end-1), 'fro') + norm(cur_p3(:, 2:end) - cur_p3(:, 1:end-1), 'fro');
for step = 1:iterations
    fprintf('Step: %d\n', step);
    cvx_begin quiet
        variable p1(n*2, T);
        variable p2(n*2, T);
        variable p3(n*2, T);
        variable u1(n, T-1); % acceleration in x stacked on acceleration in y
        variable u2(n, T-1);
        variable u3(n, T-1);
        
        p1(:, 1) == p1init;
        p2(:, 1) == p2init;
        p3(:, 1) == p3init;

        % Dynamics constraints
        for t = 1:T-1
            p1(:, t+1) == A_d * p1(:,t) + B_d * u1(:,t);
            p2(:, t+1) == A_d * p2(:,t) + B_d * u2(:,t);
            p3(:, t+1) == A_d * p3(:,t) + B_d * u3(:,t);
        end
        
        % Constraints on control
        [u1 u2 u3] <= u_max;
        [u1 u2 u3] >= u_min;
      
        
        % Objective function
        
        % collision penalty
        cost_fcn = 0;
        for t = 1:T
            dif = cur_p1(1:n, t) - cur_p2(1:n, t);
            %cost_fcn = cost_fcn + collision_scalar * max(D - ((dif' * (p1(1:n, t) - cur_p1(1:n, t)) - dif' * (p2(1:n, t) - cur_p2(1:n, t)))/norm(dif, 2) + norm(dif)), 0);
            D - ((dif' * (p1(1:n, t) - cur_p1(1:n, t)) - dif' * (p2(1:n, t) - cur_p2(1:n, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        % p2 and p3 collision
        for t = 1:T
            dif = cur_p2(1:n, t) - cur_p3(1:n, t);
            % cost_fcn = cost_fcn + collision_scalar * max(D - ((dif' * (p2(1:n, t) - cur_p2(1:n, t)) - dif' * (p3(1:n, t) - cur_p3(1:n, t)))/norm(dif, 2) + norm(dif)), 0);
            D - ((dif' * (p2(1:n, t) - cur_p2(1:n, t)) - dif' * (p3(1:n, t) - cur_p3(1:n, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        % p1 and p3 collision
        for t = 1:T
            dif = cur_p1(1:n, t) - cur_p3(1:n, t);
            % cost_fcn = cost_fcn + collision_scalar * max(D - ((dif' * (p1(1:n, t) - cur_p1(1:n, t)) - dif' * (p3(1:n, t) - cur_p3(1:n, t)))/norm(dif, 2) + norm(dif)), 0);
            D - ((dif' * (p1(1:n, t) - cur_p1(1:n, t)) - dif' * (p3(1:n, t) - cur_p3(1:n, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        
        % cost at each time step
        for t = 1:T-1
            cost_fcn = cost_fcn + (p1(:, t) - p1final)'*Q*(p1(:, t) - p1final);
            cost_fcn = cost_fcn + (p2(:, t) - p2final)'*Q*(p2(:, t) - p2final);
            cost_fcn = cost_fcn + (p3(:, t) - p3final)'*Q*(p3(:, t) - p3final);
            cost_fcn = cost_fcn + u1(:,t)'*R*u1(:,t) + u2(:,t)'*R*u2(:,t) + u3(:,t)'*R*u3(:,t);
        end
        % final cost
        cost_fcn = cost_fcn + (p1(:, T) - p1final)'*Qf*(p1(:, T) - p1final);
        cost_fcn = cost_fcn + (p2(:, T) - p2final)'*Qf*(p2(:, T) - p2final);
        cost_fcn = cost_fcn + (p3(:, T) - p3final)'*Qf*(p3(:, T) - p3final);
        
        
        % objective
        minimize(cost_fcn);
    cvx_end
    
    cur_p1 = p1;
    cur_p2 = p2;
    cur_p3 = p3;    
    costs(step+1) = cvx_optval;
end
toc

% Plot the costs
figure();
plot(costs);
title('Costs');
xlabel('Iteration');
ylabel('costs');

% Visualize the Trajectories

% % MOVIE: Replace p1, p2, p3 with your trajectories
% % ----------------------------------------------------------
% % uniform speed (colliding) trajectories
% p1 = repmat(p1init,1,T)+(p1final-p1init)*[0:1/(T-1):1]; 
% p2 = repmat(p2init,1,T)+(p2final-p2init)*[0:1/(T-1):1]; 
% p3 = repmat(p3init,1,T)+(p3final-p3init)*[0:1/(T-1):1]; 
p1 = p1(1:2, :);
p2 = p2(1:2, :);
p3 = p3(1:2, :);

[x,y,z] = cylinder(D/2,200);
x = x(1,:); y = y(1,:);
figure; cla; hold on; 
axis([-1.5,1.5,-1.5,1.5]); axis square; 
for t = 1:T
    cla;
    plot(p1(1,:),p1(2,:),'bd');
    plot(p2(1,:),p2(2,:),'rd');
    plot(p3(1,:),p3(2,:),'kd');
    plot(p1(1,t)+x,p1(2,t)+y,'b','linewidth',1.5);
    plot(p2(1,t)+x,p2(2,t)+y,'r','linewidth',1.5);
    plot(p3(1,t)+x,p3(2,t)+y,'k','linewidth',1.5);
    rectangle('Position',[-1,-1,2,2],'Linestyle','--');
    pause(0.2);
end
