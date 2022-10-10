% Constrained Model Predictive Control Function
% A: A matrix of the linear considered dynamic system
% B: B matrix of the linear considered dynamic system
% Q: status weight into the MPC cost function
% R: inputs weight into the MPC cost function
% S: final state weight (prediction horizon instant time) into the MPC cost function
% N: prediction horizon
% umin: inputs lower limit (scalar)
% umax: inputs upper limit (scalar)
% X: measured status at the current instant time
function u = mympc(A,B,Q,R,S,N,umin,umax,x)
    m=size(B,2);
    n=size(A,1);
    
    %% Q and R matrices for Open-Loop MPC (with Q1)
    Qsig = blkdiag(kron(eye(N-1),Q),S);
    Rsig = kron(eye(N),R);

    %% A and B matrices for the Open-Loop MPC
    % A matrix
    Asig = A;
    for i = 2:N
        Asig = [Asig; A^i];
    end

    % B matrix
    Bsig = [];
    temp = [];
    for i = 1:N
        temp = zeros(size(B,1)*(i-1),size(B,2));
        for j = 0:N-i
            temp = [temp; A^(j)*B];
        end
        Bsig = [Bsig temp];
    end

    %% H, F 
    H = Bsig'*Qsig*Bsig + Rsig;
    F = Asig'*Qsig*Bsig;
    ft = x'*F;
    f=ft';

    %input and status constraints definition
    lb = [repmat(umin, N*m,1)];
    ub = [repmat(umax, N*m,1)];

    options = optimset('Algorithm', 'interior-point-convex','Diagnostics','off', ...
        'Display','on');
    %solve the quadratic programming problem
    U = quadprog(H,f,[],[],[],[],lb,ub,[],options);

    %get the optimal input value (the receding horizon principle is applied)
    u = U(1:m);
end