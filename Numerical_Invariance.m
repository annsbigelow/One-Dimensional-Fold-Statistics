clear

% Set up grid
n = 80;
h = 1/n;

% Initial condition: constant function which integrates to 1
f = 1 + zeros(1,n+1);
grd=(0:h:1);

% Track the error
steps = 4;
x=grd;
f_analytical = (48.* 0.9.*(9 + 2.*x - x.^2 - 2.* x.^3 + x.^4) + ...
    1.25* (64 - 132.* x - 25.* x.^2 + 272.* x.^3 - 31.* x.^4 - 126.* x.^5 + ...
    42.* x.^6))./(96.*(-2 + x).^2 .*(1 + x).^2);
errvec = zeros(1,steps);

% Iterate f
%for m = 1:15
    %Nvec(m)=n; % Increase n to check convergence
    for i = 1:steps
        f = fNew(f,n,grd);
        %disp(['Integral value:', num2str(trapz(grd,f))]);
        %errvec(i)=sqrt(h*sum((f_analytical-f).^2)); % Weighted norm
        %errvec(i) = norm(f_analytical-f);
    end
%     errvec(m)=sqrt(h*sum((f_analytical-f).^2)); % Weighted norm
%     n=n*2;
%     h=1/n;
%     f = 1 + zeros(1,n+1);
%     grd=(0:h:1);
%     x=grd;
%     f_analytical = (48.* 0.9.*(9 + 2.*x - x.^2 - 2.* x.^3 + x.^4) + 1.25* (64 - 132.* x - 25.* x.^2 + 272.* x.^3 - 31.* x.^4 - 126.* x.^5 + 42.* x.^6))./(96.*(-2 + x).^2 .*(1 + x).^2);
%end

%Plot to check the final iteration
plot(grd,f);
xlabel('s');
ylabel('f(s)');
hold on; 
plot(grd, f_analytical);
legend('numerical','analytical');
hold off; 


%[C,p]=my_loglog(Nvec,errvec,'n','Error','Error for Increasing Points');

function [const,expo] = my_loglog(xin,yin,my_xlabel,my_ylabel,my_title)
    % Fit a linear to error
    pol = polyfit(log(xin), log(yin), 1);
    const = exp(pol(2));
    expo = pol(1);
    
    % Generate y values for the fitted line
    yfit = polyval(pol, log(xin));
    
    % Plot 
    loglog(xin, yin, 'o');
    hold on;
    %loglog(xin, exp(yfit), '-');
    xlabel(my_xlabel);
    ylabel(my_ylabel);
    title(my_title);
    hold off;
end
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
    I1 = integral(@(t) t.*L(s.*t,f,n,grd), 1/(2-s), 1);
    I2 = integral(@(t) t.*(L(s.*t,f,n,grd) + L(t.*(2-s),f,n,grd)), 1/2, 1/(2-s));
    z = I1 + I2;
end

function fnew = fNew(f,n,grd)
    fnew = zeros(1,n+1);
    for i = 1:n+1
        s = grd(i);
        fnew(i) = F(s,f,n,grd) + F(1-s,f,n,grd);
    end
end