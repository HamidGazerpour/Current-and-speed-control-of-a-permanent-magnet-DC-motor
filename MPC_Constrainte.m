function [ u ] = MPC_Constrainte( A,B, Q, R, S, N, umin, umax, xmin, xmax, x)

options = optimset('quadprog');
options = optimset(options, 'LargeScale', 'off', 'Display' , 'off');
rho = 100;

Ac = [];
for i = 1:N
    Ac = [Ac; A^i];
end

Bc = [];
for i= 1:N
    line = [];
    for j = 1:N
        if(j>i)
            el = zeros(size(B));
        else
            el = A^(i-j) * B;
        end
        line = [line el];
    end
    Bc = [Bc;line];
end

Qc = kron(eye(N-1), Q);
Qc = blkdiag(Qc, S);

Rc = kron(eye(N), R);
H = 2 * (Bc' * Qc * Bc + Rc);
H = blkdiag(H, 0);
f = 2 * x' * Ac' * Qc * Bc;
f = [f rho];

umin_constr = [umin .* ones(N, 1); 0];
umax_constr = [umax .* ones(N, 1); inf];

[A_state, b_state] = build_state_contraint_matrix(Ac, Bc, N, xmin, xmax, x);

A_leq = [A_state];
b_leq = [b_state];


[u, cost, flag] = quadprog(H, f, A_leq, b_leq, [], [], umin_constr, umax_constr, [], options);
handle_error(flag);
cost = cost - u(11) * rho + x' * Ac' * Qc * Ac * x;
u = [u(1)];
end

function handle_error(flag)
if(flag == -3)
    error('Unbounded optimisation problem. Execution aborted.');
end
if(flag == -2)
    error('Infeasible optimisation problem. Execution aborted.');
end
if(flag == 0)
    error('Maximum number of iteration exceeded. Execution aborted.')
end
end

function [A, b] = build_state_contraint_matrix(Ac, Bc, N, x_min, x_max, x)
ev = repmat([0 1]', N, 1);
x_max_vect = repmat(x_max, N, 1);
b = x_max_vect - Ac * x;
A = [Bc -ev];

x_min_vect = repmat(x_min, N, 1);
A = [A; -[Bc ev]];
b = [b; - x_min_vect + Ac * x];
end
