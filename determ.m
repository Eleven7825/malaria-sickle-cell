% compatibility: Matlab R 2019b+

function [lambda,S,I,R,N] = determ(ops)

% argument default and validations
arguments
    ops.N = 1e5;
    ops.I = 2e4;
    ops.R = 3e4;
    ops.T = 600;
    ops.d =0.01
    ops.d_ = 0.02
    ops.b  = 0.02
    ops.b_ = 0.01
    ops.a = 0.5;
    ops.Plot = 1;
    ops.saveVariable = 1
end

b = ops.b; d = ops.d; b_ = ops.b_; d_ = ops.d_;
N = ops.N; I = ops.I; R = ops.R; S = N-I-R;
a = ops.a; T = ops.T;

if ~(b>b_)
    warning("b must be greater than b_")
elseif ~(d<d_)
    warning("d must be less than d_")
elseif ~(b>d)
    warning("b must be greater than d")
elseif ~(b_<d_)
    warning("b_ must be less than d_")
end

if ops.saveVariable==1
    Ssave = zeros(1,T+1); Isave = zeros(1,T+1); Rsave = zeros(1,T+1);Nsave = zeros(1,T+1);
    lambdas = zeros(1,T+1);
    Ssave(1) = S; Isave(1) = I; Rsave(1) = R; Nsave(1) = N;
end

% the main loop for the simualtion 
for t = 1: T
    pgtb00 = (b*S+b_*I); % the 00 type people going to born (pgtb00)
    pgtb01 = (b*R);         % the 01 type people going to born (pgtb01)
    P00_00 = (S+I)/N+R/(2*N); % probability for 00 type woman give birth to 00 type child
    P00_01 = 1-P00_00 ;% probability for 00 type woman give birth to 01 type child
    P01_00 = (S+I)/(2*N)+R/(4*N); % probability for 01 type woman give birth to 00 type child
    P01_01 = (S+I)/(2*N)+R/(2*N); % probability for 01 type woman give birth to 01 type child
    BS = pgtb00*P00_00+pgtb01*P01_00;% new suspetible baby(00)
    BR = pgtb00*P00_01+pgtb01*P01_01;% new resistant baby(01)
    
    % update S I R 
    % we are suposing dt=1 so I didn't create a variable called dt
    S = S+BS-d*S-a*I*S/N;
    I = I-d_*I+a*I*S/N;
    R = R+BR-d*R;
    
    % if any quatity is smaller than 0, make it equal to 0
    if S<0; S=0; end 
    if I<0; I=0; end
    if R<0; R = 0; end
    
    N = S+I+R;
    
    % save the variable
    if ops.saveVariable == 1
        Ssave(t+1) = S; Isave(t+1) = I; Rsave(t+1) = R;
        
    end
    Nsave(t+1) = N;
    lambdas(t+1) = (N-Nsave(t))/N;
end

lambda = mean(lambdas(floor(ops.T/2):end));

% plot SIR if speficied
if ops.Plot == 1
    hold on
    plot(1:T+1,Ssave,'LineWidth',1.5);plot(1:T+1,Isave,'LineWidth',1.5);plot(1:T+1,Rsave,'LineWidth',1.5);
    set(gca,'YScale','log')
    legend('S','I','R')
    xlabel('Period')
    ylabel('Number')
    title(sprintf('Deterministic model with a = %2f',a));
    hold off
end

% save the variables if required
if ops.saveVariable==1
    save('SIR','Ssave','Isave','Rsave','Nsave','lambdas')
end
end