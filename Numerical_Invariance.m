clear

% Set up grid
n = 50;
h = 1/n;

% Initial condition: constant function which integrates to 1
f = 1 + zeros(1,n+1);

grd=(0:h:1);

% Iterate f
steps = 10;
for i = 1:steps
    f = fNew(f,n,grd);
    f = f./trapz(f); % This normalization step shouldn't be included in our final version because it is cheating.
    % Comment/uncomment line 16 to see the true integral value of f. Currently (without line 16),
    %   f is not properly normalized.
    disp(['Integral value:', num2str(trapz(f))]);
end

% Plot to check the final iteration
plot(grd,f);

xlabel('s');
ylabel('f(s)');

% Functions defined below.

function z = L(x,f,n,grd) % Lagrange polynomial interpolant for each integration step
    % X,Y are sized k+1
    % Choose k points from 0 to 1 (not perfectly equally-spaced)
    k = 7; % Max degree
    space = floor(n/k);
    j = 1;
    X = zeros(1,k+1);
    Y = zeros(1,k+1);
    for i = 1:k+1
        X(i) = grd(j);
        Y(i) = f(j); 
        j = j + space;
    end 

    % Lines 41-43 yield Runge Phenom. (or choosing any more points than ~7)
%     k=n;
%     Y=f;
%     X=grd;

    z = 0;
    for j = 1:k+1
        prod = 1; 
        for m = 1:k+1
            if m ~= j
                prod = prod.*(x-X(m))./(X(j)-X(m)); % x is real
            end
        end
        z = z + Y(j).*prod;
    end
end


function z = F(s,f,n,grd) 
    % Trapezoid method for integration using Lagrange polynomials (this
    %   expression was found via paper & pen)
    % z = (s/4)*((1/2)*L(s/2,f,n,grd) + ((s+1)/2)*L(s*(s+1)/2,f,n,grd))+((1-s)/4)*(((s+1)/2)*(L(s*(s+1)/2,f,n,grd)+L((s+1)*(2-s)/2,f,n,grd)) + L(s,f,n,grd) + L(2-s,f,n,grd));
    
    % Try MATLAB's "integral()". Works!
    I1 = integral(@(t) t.*L(s.*t,f,n,grd), 1/2, (s+1)/2);
    I2 = integral(@(t) t.*(L(s.*t,f,n,grd) + L(t.*(2-s),f,n,grd)), (s+1)/2, 1);
    z = I1 + I2;

    % Try MATLAB's built-in trapz... no good
%     t1=linspace(0.5,(s+1)/2,25);
%     t2=linspace((s+1)/2,1,25);
%     for i = 1:25
%         Y1(i)=t1(i)*L(s*t1(i),f,n,grd);
%         Y2(i)=t2(i)*(L(s*t2(i),f,n,grd)+L(t2(i)*(2-s),f,n,grd));
%     end
%     Y=Y1+Y2;
%     z=trapz(Y);
end

% function z = F(s,f,grd) % Trap. method using linear interpolation
%     % xq are "query" points where we want the interpolant to be evaluated.
%     xq = [s/2, s*(s+1)/2, (s+1)*(2-s)/2, s, 2-s];
%     grd(end) = 2; % Extend last point all the way out to x=2. Seems sketchy.
%     p = interp1(grd,f,xq); % Interpolate f at ALL the grid points and evaluate at xq
% 
%     z = (s/4)*((1/2)*p(1) +((s+1)/2)*p(2))+((1-s)/4)*(((s+1)/2)*(p(2)+p(3)) + p(4) + p(5));
% end

function fnew = fNew(f,n,grd)
    fnew = zeros(1,n+1);
    for i = 1:n+1
        s = grd(i);
        % Arguments of function "F()" have to be changed according to
        %    linear/Lagrange interpolation
        fnew(i) = F(s,f,n,grd) + F(1-s,f,n,grd);
    end
end