% compatibility: Matlab R 2019b+

function [lambda,S,I,R,N] = nondeterm(ops)

% argument default and validations

arguments
    ops.N = 1e5;
    ops.I = 2e4;
    ops.R = 3e4;
    ops.T = 600;
    ops.d=0.01
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
     pgtb00 = (b*S+b_*I); % the 00 type woman that are going to born
    pgtb01 = (b*R);         % the 01 type woman that are going to born
    P00_00 = (S+I)/N+R/(2*N); % probability for 00 type woman give birth to 00 type child
    P00_01 = 1-P00_00 ;% probability for 00 type woman give birth to 01 type child
    P01_00 = (S+I)/(2*N)+R/(4*N); % probability for 01 type woman give birth to 00 type child
    P01_01 = (S+I)/(2*N)+R/(2*N); % probability for 01 type woman give birth to 01 type child
    BS = pgtb00*P00_00+pgtb01*P01_00;% new suspetible baby(00)
    BR = pgtb00*P00_01+pgtb01*P01_01;% new resistant baby(01)
    
    % dealing with infected people
    % number of Infected new borns
    Ib = sum(rand([1 I],'single')<b_);%number of Infected people get borned
    
    % Resistant baby born from Infected mom
    I2R =  sum(rand([1 Ib], 'single')<P00_01);
    
    % Suspectible baby born from Infected mom
    I2S = sum(rand([1 Ib], 'single')<P00_00);
    
    % dealing with suspectible people
    Sb = sum(rand([1 S], 'single')<b);
    
    % Resistant baby born from suspectible mom
    S2R =  sum(rand([1 Sb], 'single')<P00_01); 
    
    % Suspectible baby born from suspectible mom
    S2S =  sum(rand([1 Sb], 'single')<P00_00); 
    
    % dealing with Resistant people
    Rb = sum(rand([1 R],'single')<b);
    
    % Resistant baby born from Resistant mom
    R2R =  sum(rand([1 Rb], 'single')<P01_01); 
    
    % Suspectible baby born from resistant mom
    R2S =  sum(rand([1 Rb], 'single')<P01_00); 
    
    % consoder the spread of the desease
    S2I = sum(rand([1 S], 'single')<a*I/N);
    
    % consider some people dies
    Id = sum(rand([1 I], 'single')<d_);
    Sd = sum(rand([1 S], 'single')<d);
    Rd = sum(rand([1 R], 'single')<d);
    
    S = S + I2S + R2S + S2S - Sd-S2I;
    I = I+S2I-Id;
    R = R +I2R + R2R +S2R - Rd;
    
    % if any quatity is smaller than 0, make it equal to 0
    if S<0; lambda = 0;return; end
    if I<0; lambda = 0;return; end
    if R<0; lambda = 0;return; end
    
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
    plot(1:T+1,Ssave);plot(1:T+1,Isave);plot(1:T+1,Rsave);
    set(gca,'YScale','log')   
    xlabel('Period')
    ylabel('Number')
    title(sprintf('Non deterministic model with a = %2f',a));
    hold off
end

% save the variables if required
if ops.saveVariable==1
    save('SIR','Ssave','Isave','Rsave','Nsave','lambdas')
end