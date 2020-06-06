function [p1_ret, p2_ret, p3_ret, costs] = traj_opt_fcn(p1init, p2init, p3init, p1final, p2final, p3final, cur_p1, cur_p2, cur_p3, iterations)

n = 2; N = 3; T = 30; D = 0.4;

costs = zeros(1, iterations+1);
costs(1) = norm(cur_p1(:, 2:end) - cur_p1(:, 1:end-1), 'fro') + norm(cur_p2(:, 2:end) - cur_p2(:, 1:end-1), 'fro') + norm(cur_p3(:, 2:end) - cur_p3(:, 1:end-1), 'fro');
for step = 1:iterations
    fprintf('Step: %d\n', step);
    cvx_begin quiet
        variable p1(n, T);
        variable p2(n, T);
        variable p3(n, T);
        
        p1(:, 1) == p1init;
        p2(:, 1) == p2init;
        p3(:, 1) == p3init;
        
        p1(:, end) == p1final;
        p2(:, end) == p2final;
        p3(:, end) == p3final;

        % p1 and p2 collision
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
        minimize(norm(p1(:, 2:end) - p1(:, 1:end-1), 'fro') + norm(p2(:, 2:end) - p2(:, 1:end-1), 'fro') + norm(p3(:, 2:end) - p3(:, 1:end-1), 'fro'));
    cvx_end
    
    cur_p1 = p1;
    cur_p2 = p2;
    cur_p3 = p3;    
    costs(step+1) = norm(cur_p1(:, 2:end) - cur_p1(:, 1:end-1), 'fro') + norm(cur_p2(:, 2:end) - cur_p2(:, 1:end-1), 'fro') + norm(cur_p3(:, 2:end) - cur_p3(:, 1:end-1), 'fro');
end

p1_ret = p1;
p2_ret = p2;
p3_ret = p3;
