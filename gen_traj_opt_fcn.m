function [P_ret] = gen_traj_opt_fcn(Pinit, Pfinal, cur_P, iterations)
shape = size(cur_P);
n = shape(1); N = shape(3); T = shape(2); D = 0.4;

costs = zeros(1, iterations+1);
for i = 1:N
    costs(1) = costs(1) + norm(cur_P(:, 2:end, i) - cur_P(:, 1:end-1,i), 'fro')
end

for step = 1:iterations
    fprintf('Step: %d\n', step);
    cvx_begin quiet
        variables P(n,T,N) p1(n,T) p2(n,T) p3(n,T);  
        
        % initial and final constraint
        reshape(P(:,1,:),n,N) == Pinit; 
        reshape(P(:,end,:),n,N) == Pfinal;         

        % collision constraint
%         for i = 1:N
%             for j = 1:N
%                 if i > j
%                     for t = 1:T
%                         dif = cur_P(:, t, i) - cur_P(:, t, j);
%                         D - ((dif' * (P(:, t, i) - cur_P(:, t, i)) - dif' * (P(:, t, j) - cur_P(:, t, j)))/norm(dif, 2) + norm(dif)) <= 0;
%                     end
%                 end
%             end
%         end        

% p1 and p2 collision
        p1 = P(:,:,1); p2 = P(:,:,2); p3 = P(:,:,3);
        cur_p1 = cur_P(:,:,1);
        cur_p2 = cur_P(:,:,2);
        cur_p3 = cur_P(:,:,3);
        for t = 1:T
            dif = cur_p1(:, t) - cur_p2(:, t);
            D - ((dif' * (p1(:, t) - cur_p1(:, t)) - dif' * (p2(:, t) - cur_p2(:, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        % p2 and p3 collision
        for t = 1:T
            dif = cur_p2(:, t) - cur_p3(:, t);
            D - ((dif' * (p2(:, t) - cur_p2(:, t)) - dif' * (p3(:, t) - cur_p3(:, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        % p1 and p3 collision
        for t = 1:T
            dif = cur_p1(:, t) - cur_p3(:, t);
            D - ((dif' * (p1(:, t) - cur_p1(:, t)) - dif' * (p3(:, t) - cur_p3(:, t)))/norm(dif, 2) + norm(dif)) <= 0;
        end
        
        % objective
        J = 0;
        for i = 1:N
            J = J + norm(cur_P(:, 2:end, i) - cur_P(:, 1:end-1, i), 'fro');
        end
        minimize (J)
    cvx_end
    
    % update cur_P 
    cur_P = P
    
    % update cost
    for i = 1:N
        costs(step+1) = costs(step+1) + norm(cur_P(:, 2:end, i) - cur_P(:, 1:end-1,i), 'fro');
    end
end

P_ret = P;
