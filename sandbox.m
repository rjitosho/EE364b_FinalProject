close all;
clear all;
tic

% Problem Definition
n = 2; N = 20; T = 30; D = 0.3; T_s = 0.1;

% Define our cost matrices
Q = diag([.1 .1 .1 .1]) * 1.0; % cost associated with the state
Qf = eye(2*n,2*n) * 100.0; % cost associated with the final state
R = eye(n, n) * 0.1; % cost associated with the control effort 
collision_scalar = 10.0;

% Dynamics matrices
A_d = [eye(n, n) T_s*eye(n,n); zeros(n,n) eye(n,n)];
B_d = [T_s^2/2*eye(n,n); T_s*eye(n,n)];

% Control constraints
u_max = 20.0;
u_min = -20.0;

% State constraints
v_min = -D/(2*T_s);
v_max = D/(2*T_s);

% Initial and final positions
Pinit = zeros(2*n, N);
radius = 1;
for i = 1:N
    theta = (i-1)/N * 2 * pi;
    cur_x = radius * cos(theta);
    cur_y = radius * sin(theta);
    Pinit(:, i) = [cur_x; cur_y; 0; 0];
end

locations = 1:N;
end_locations = randperm(N);
while(all(locations == end_locations))
    end_locations = randperm(N);
end
Pfinal = Pinit(:, end_locations);

% Initial trajectory
cur_P = zeros(2*n, T, N);
for i = 1:N
    cur_P(:, :, i) = repmat(Pinit(:, i), 1, T) + (Pfinal(:, i) - Pinit(:, i)) * [0:1/(T-1):1];
end

%% SCP 
% % OLD INITIAL TRAJ
% p1init = [-1,0, 0, 0]'; p2init = [-1,0.5, 0, 0]'; p3init = [-1,1, 0, 0]';
% p1final = [1,0, 0, 0]'; p2final = [1,-0.5, 0, 0]'; p3final = [1,-1, 0, 0]';
% cur_P = zeros(2*n,T,N);
% cur_P(:,:,1) = repmat(p1init,1,T)+(p1final-p1init)*[0:1/(T-1):1];
% cur_P(:,:,2) = repmat(p2init,1,T)+(p2final-p2init)*[0:1/(T-1):1]; 
% cur_P(:,:,3) = repmat(p3init,1,T)+(p3final-p3init)*[0:1/(T-1):1]; 
% Pinit = [p1init p2init p3init];
% Pfinal = [p1final p2final p3final];

% calculate the velocities for initial traj
for i = 1:N
    cur_P(n+1:2*n,:,i) = [(cur_P(1:n,2:T,i)-cur_P(1:n,1:T-1,i))/T_s [0;0]];   
end

iterations = 10;

costs = zeros(1, iterations+1);
for i = 1:N
    costs(1) = costs(1) + norm(cur_P(:, 2:end, i) - cur_P(:, 1:end-1,i), 'fro');
end

for step = 1:iterations    
    fprintf('Step: %d\n', step);
    cvx_begin quiet
        variables P(2*n,T,N) U(n,T,N);     
        
        % initial constraint
        reshape(P(:,1,:),2*n,N) == Pinit; 
%         reshape(P(:,end,:),2*n,N) == Pfinal;
        
        % Dynamics constraints
        for t = 1:T-1
            for i = 1:N
                P(:, t+1, i) == A_d * P(:,t,i) + B_d * U(:,t,i);
            end
        end      
        
        % collision constraint
        for i = 1:N
            for j = i+1:N               
                for t = 1:T
                    dif = cur_P(1:n, t, i) - cur_P(1:n, t, j);
                    D - ((dif' * (P(1:n, t, i) - cur_P(1:n, t, i)) - dif' * (P(1:n, t, j) - cur_P(1:n, t, j)))/norm(dif, 2) + norm(dif)) <= 0;
                end                
            end
        end  
        
        cost_fcn = 0;        
        % cost at each time step
        for i = 1:N
            for t = 1:T-1
                cost_fcn = cost_fcn + (P(:, t, i) - Pfinal(:,i))'*Q*(P(:, t, i) - Pfinal(:,i));
                cost_fcn = cost_fcn + U(:,t,i)'*R*U(:,t,i);
            end
            % final cost
            cost_fcn = cost_fcn + (P(:, T, i) - Pfinal(:,i))'*Qf*(P(:, T, i) - Pfinal(:,i));
        end
        
        % objective
        minimize(cost_fcn);
    cvx_end
    
    % update cur_P 
    cur_P = P;
    
    % update cost
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
[x,y,z] = cylinder(D/2,200);
x = x(1,:); y = y(1,:);
figure; cla; hold on; 
axis([-1.5,1.5,-1.5,1.5]); axis square; 
vehicle_colors = rand(3,N);
for k = 1:10
    for t = 1:T
        cla;    
        for i = 1:N
            plot(P(1,:,i),P(2,:,i),'d','color',vehicle_colors(:,i));
            plot(P(1,t,i)+x,P(2,t,i)+y,'color',vehicle_colors(:,i),'linewidth',1.5);
        end

        rectangle('Position',[-1,-1,2,2],'Linestyle','--');
        pause(0.2);
    end
end
