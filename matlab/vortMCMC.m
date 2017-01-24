% Implementation of the vortical MCMC follows.

function x = vortMCMC(Q,q,p,x0,a,N)
%{
        Markov Chain Monte Carlo with vortices.
        
        The following code uses the proposal transition matrix Q to construct a non-reversible Markov chain
        X_n with invariant distribution p. To do this, a modification of the classical acceptance probability is
        proposed in [1]. This modification adds an extra part called the vorticity matrix G. In this version of the code,
        the vorticity matrix used is,
    
            G(x,y) = q(x,y)Q(x,y)-q(y,x)Q(y,x)
            Ga(x,y) = a*G(x,y)
        
        Input:
             Q : (nxn matrix) : A transition probability matrix. Used to make proposals.
             q : (nx1 vector) : Stationary distribution corresponding to Q.
             p : (nx1 vector) : Target distribution to sample from (unnormalized).
             x0: (1x1 integer): Initial point to start from.
             a : (1x1 float)  : Parameter controling the vorticity matrix.
             N : (1x1 integer): Number of samples to create.
                     
        Output:
            x: N samples from the invariant distribution.
        
        Warning :
            Plots need to be made to make sure that x has equilibrated to p.
            Then, an appropriate number of samples can be thrown away.
    
%}

M = size(Q,1); % Dimension of Q. Here Q is assumed to be a square matrix.

% Construction of vorticity matrix associated with Q and q.
G = zeros(M,M);
for i = 1:M
    for j = 1:M
        G(i,j) = a*(q(i)*Q(i,j)-q(j)*Q(j,i)); % Matrix elements are scaled by a.
    end
end


x = zeros(1,N); % Array to store samples
x(1) = x0;

for i = 1: N
    % Generate proposal according to Q
    probs = cumsum(Q(x(i),:)); % Probabilities for this step.
    u = rand();
    for j = 1:M
        if(u<probs(j))
            prop = j;
            break;
        end
    end

    % Calculate acceptance ratio
    acc_ratio = (G(x(i),prop)+p(prop)*Q(prop,x(i)))/(p(x(i))*Q(x(i),prop));

    % Acceptance probability
    acc_prob  = min(1,acc_ratio);

    u = rand();

    if(u<acc_prob)% Accept sample.
        x(i+1) = prop;
    else % Reject sample.
        x(i+1) = x(i);
    end
    
end

   
end