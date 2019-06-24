%----------------------------------------------------------------------------------------------------
% Criteria implementing optimality conditions for unconstrained
% optimization problems
%----------------------------------------------------------------------------------------------------

clear all;
close all;

N=2; % number of variables (edit this according to objective function)
for i=1:N 
    syms("x"+i);
end
    
% objective function
f = 8*x1^2 + 3*x1*x2 + 7*x2^2 - 25*x1 + 31*x2 - 29
%f= x1^2 + x1*x2 + 1.5*x2^2 - 2*log(x1) - log(x2)
%f= x1^2 - 2*x1*x2^2 + x2^4 - x2^5
%f = x1^2 + 2*x2^2 + 5*x3^2 - 2*x1*x2 - 4*x2*x3 - 2*x3

% gradient of the objective function
for i=1:N
      GradF(i,1)=diff(f, "x"+i);
end
disp("Gradient of the f= ");
disp(GradF);

% Hessian matrix of the objective function
for i=1:N
    for j=1:N
      HF(i,j)=diff(diff(f,"x"+j), "x"+i); 
    end
end
disp("Hessian Matrix= ");
disp(HF);                                                                              

% Solution of the system gradient f=0  and determination of the stationary points of the objective function f
xi=sym(zeros(1,N));
xii=sym(zeros(1,N));
for i=1:N
     GradF2(1,i) = GradF(i,1)==0;
     xi(1,i) = ("x"+i);
     xii(i,1) = ("x"+i);
end

X = solve(GradF2, xi);
x_star_num = size(X.x1,1); % get number of solutions

% Displays all x stars
for j=1:x_star_num
    for i=1:N
         x_star(j,i) = X.("x"+i)(j) ;
    end
    disp("x_star_"+j+"=");
    disp(x_star(j,:));
end

for j=1:x_star_num
    HF_star=subs(HF,xi,x_star(j,:));   % Hessian Matrix of a certain x star
    disp("Hessian Matrix of x_star_"+j+"=");
    disp(HF_star);
    
    % Computes & displays Quadratic form of a certain x star
    Qf = xi * HF_star * xii;
    disp("******************************");
    disp("* Quadratic form of x_star_"+j+" *");
    disp("******************************");
    disp(Qf(1,1));
   
    % Computes & displays Determinants of a certain x star
    D_star(1,1) = HF_star(1,1);
    if N==2
        D_star(2,1)=det(HF_star);
    end
    if N>2
        for i=2:N
            array = HF_star(1:i,1:i);
            D_star(i,1) = det(array);
        end
            
    end
    disp("****************************");
    disp("* Determinants of x_star_"+j+" *");
    disp("****************************");
    pos=0; neg=0; zero=0;
    for i=1:N
         fprintf('Δ%d = %d \n',i , D_star(i,1));  
         if int8(D_star(i,1)) > 0 % If Δκ is positive
             pos = pos + 1;
         end
         if rem(i, 2) == 0 % If k is even number
             if int8(D_star(i,1)) > 0 % If Δκ is positive
                 neg = neg + 1;
             end
         else               % If k is odd number
             if int8(D_star(i,1)) < 0 % If Δκ is negative
                 neg = neg + 1;
             end
         end
         if int8(D_star(i,1)) == 0 % If Δκ is zero
             zero = 1;
         end
    end
    if zero == 1
         disp("There is at least 1 Determinant = 0 -> No conclusion for x_star_"+j);
    elseif pos == N
        disp("All Determinants are positive -> x_star_"+j+ " is local minimizer");
    elseif neg == N
        disp("Determinants Δκ are negative for κ=odd number & Δκ are positive for κ=even number -> x_star_"+j+ " is local maximizer");
    else
        disp("At least one Determinant Δκ is negative for κ=even number -> x_star_"+j+ " is not local optimizer");
    end
    fprintf('\n');

    
    % Computes & displays Eigenvalues of a certain x star
    eig_star = eig(HF_star);
    disp("***************************");
    disp("* Eigenvalues of x_star_"+j+" *");
    disp("***************************");
    pos=0; neg=0; zero=0;
    for i=1:N
         fprintf('λ%d = %d \n',i , eig_star(i,1));
         if int8(eig_star(i,1)) > 0       % If λκ is positive
             pos = pos + 1;
         elseif int8(eig_star(i,1)) < 0   % If λκ is negative
             neg = neg + 1;
         else                       % If λκ is zero
             zero = 1;
         end
    end
    if zero == 1
        disp("There is at least 1 Eigenvalue = 0 -> No conclusion for x_star_"+j);
    elseif pos == N
        disp("All Eigenvalues are positive -> x_star_"+j+ " is local minimizer");
    elseif neg == N
        disp("All Eigenvalues are negative -> x_star_"+j+ " is local maximizer");
    else
        disp("There are at least 1 posite & 1 negative Eigenvalue -> x_star_"+j+ " is not local optimizer");
    end
    fprintf('\n');
end

if N == 2
    ezsurf(f,[-0.01  0.01]);
end

