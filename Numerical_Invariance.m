clear

% Set up grid
n = 50;
h = 1/n;

% Initial condition: constant function which integrates to 1
f = 1 + zeros(1,n+1);
grd=(0:h:1);
disp(['Integral value:', num2str(trapz(grd,f))]);

% Iterate f
steps = 10;
for i = 1:steps
    f = fNew(f,n,grd);
    %f = f./trapz(f); % This normalization step shouldn't be included in our final version because it is cheating.
    % Comment/uncomment line 16 to see the true integral value of f. 
    disp(['Integral value:', num2str(trapz(grd,f))]);
end

% Plot to check the final iteration
plot(grd,f);
xlabel('s');
ylabel('f(s)');


function z = L(x,f,n,grd) % Lagrange polynomial interpolant for each integration step
    % X,Y are sized k+1
    % Choose k points from 0 to 1 (not perfectly equally-spaced)
    k = 4; % Max degree
    space = floor(n/k);
    j = 1;
    X = zeros(1,k+1);
    Y = zeros(1,k+1);
    for i = 1:k+1
        X(i) = grd(j);
        Y(i) = f(j); 
        j = j + space;
    end 

    % Lines 41-43 yield Runge Phenomenon (or choosing any more points than ~7)
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
    %   expression was found via paper & pen). 
    z = (1/2)*((1-s)/(2-s)*(L(s/(2-s),f,n,grd)/(2-s) + L(s,f,n,grd)) + s/(2*(2-s))*((L(s/2,f,n,grd)+L((2-s)/2,f,n,grd))/2 + (L(s/(2-s),f,n,grd)+f(end))/(2-s)));
    
    % MATLAB's "integral()"
%     I1 = integral(@(t) t.*L(s.*t,f,n,grd), 1/(2-s), 1);
%     I2 = integral(@(t) t.*(L(s.*t,f,n,grd) + L(t.*(2-s),f,n,grd)), 1/2, 1/(2-s));
%     z = I1 + I2;
end

function fnew = fNew(f,n,grd)
    fnew = zeros(1,n+1);
    for i = 1:n+1
        s = grd(i);
        fnew(i) = F(s,f,n,grd) + F(1-s,f,n,grd);
    end
end
