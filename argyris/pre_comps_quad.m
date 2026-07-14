clear
% EXAMPLE vertex of the reference triangle (0, 1, or 2)
vertex = 2;

% Compute the basis of the Argyris element on the reference triangle
% M(i,:) holds the coefficients of the i-th ref. basis function in the
%   monomial basis
M = reference();

% Plugging in (0,0), (1,0), or (0,1) will produce e_1, e_2, or
%   e_3,.. respectively for hat_z, grads
% NOTE: THERE ARE ONLY THREE OPTIONS FOR (X,Y) HERE.
% Function values and derivatives are useful for quadrature later.
if vertex == 0
    eval = monomials(0,0);
    mx = derx(0,0); my = dery(0,0);
    mxx = derxx(0,0); mxy = derxy(0,0); myy = deryy(0,0);
elseif vertex == 1
    eval = monomials(1,0);
    mx = derx(1,0); my = dery(1,0);
    mxx = derxx(1,0); mxy = derxy(1,0); myy = deryy(1,0);
else 
    eval = monomials(0,1);
    mx = derx(0,1); my = dery(0,1);
    mxx = derxx(0,1); mxy = derxy(0,1); myy = deryy(0,1);
end

% Each row/element of hat_z is phi_i evaluated at (x,y) 
hat_z = M*eval;

% First derivatives
grads = zeros(21,2);
grads(:) = [mx(:) my(:)];
% 1st col of grads is d_x phi_i evaluated at (x,y); 2nd col d_y...
grads = M*grads;

% Second derivatives
hess = zeros(21,3);
hess(:) = [mxx(:) mxy(:) myy(:)];
hess = M*hess;

%% Try n-point Gauss-Legendre quadrature
% Compute the quadrature weights
n=6;
[w,xi] = get_weights(n);

s = @(xi,b) (b/2)*xi + b/2;
%% Approximate the mass integrals
phi_eval = @(x,y) M*monomials(x,y); % Function values
S = zeros(21,21);
for a=1:21
    for b=1:21
    for j=1:n
        prefac = w(j)*(1-s(xi(j),1))/2;
        sumi=0;
            for i=1:n
                x = s(xi(j),1);
                phi = phi_eval(x, s(xi(i),1-x));
                sumi = sumi + w(i)*phi(a)*phi(b);
            end
        S(a,b) = S(a,b) + prefac*sumi;
    end
    S(a,b) = S(a,b)/2;
    end
end

S_exact = readmatrix("S_precomp.csv");
format long e;
for i=1:21
    for j=1:21
        diff = abs(S_exact(i,j)-S(i,j));
        if (diff > 10^-8)
            disp([i,j,diff]);
        end
    end
end
disp(cond(S));
%% Approximate the stiffness integrals
% BELOW IS A TEST
mxx = @(x,y) derxx(x,y); mxy = @(x,y) derxy(x,y); myy = @(x,y) deryy(x,y);

x = s(xi(1),1);
y = s(xi(1),1-x);
xx=mxx(x,y); xy=mxy(x,y); yy=myy(x,y);
hess(:) = [xx(:) xy(:) yy(:)];
hess_eval = M*hess; % Hessian

%% Output M and S 
S_out = reshape(permute(S,[2,1]),[],1);
M_out = reshape(permute(M,[2,1]),[],1);

cID = fopen('pre_comps.hh','w');
fprintf(cID,'extern const double M[%d];\n',21*21);
fprintf(cID,'extern const double S[%d];',21*21);
fclose(cID);

cID = fopen('pre_comps.cc','w');
fprintf(cID, '#include "pre_comps.hh"\n\n');
fprintf(cID, 'const double S[%d] = {\n',21*21);
for k = 1:21*21
    fprintf(cID, '%.17g', S_out(k)); % Use full double precision
    if k < 21*21
        fprintf(cID, ',');
    end
    if mod(k,8)==0
        fprintf(cID, '\n');
    end
end
fprintf(cID, '};\n');

fprintf(cID, 'const double M[%d] = {\n',21*21);
for k = 1:21*21
    fprintf(cID, '%.10g', M_out(k)); 
    if k < 21*21
        fprintf(cID, ',');
    end
    if mod(k,21)==0
        fprintf(cID, '\n');
    end
end
fprintf(cID, '};');
fclose(cID);
%%
function [w,xi] = get_weights(n)
    % Legendre polynomials
    syms x; 
    P_n = legendreP(n,x);

    % Find zeros of Legendre polynomial
    r = vpasolve(P_n==0,x,1/2);
    xi = double(r);
    xi  = unique(xi);

    dP = diff(P_n,x);
    w = zeros(n,1);
    for i=1:n
        d = subs(dP,x,xi(i));
        w(i) = 2/((1-xi(i)^2)*d^2);
    end
end

function mxx = derxx(x,y)
    mxx = [0 0 0 0 0 0 ...
            0 0 0 0 0 ...
            2 2*y 2*y^2 2*y^3 ...
            6*x 6*x*y 6*x*y^2 ...
            12*x^2 12*x^2*y 20*x^3]';
end

function mxy = derxy(x,y)
    mxy = [0 0 0 0 0 0 0 ...
            1 2*y 3*y^2 4*y^3 0 ...
            2*x 4*x*y 6*x*y^2 0 ...
            3*x^2 6*x^2*y ...
            0 4*x^3 0]';
end

function myy = deryy(x,y)
    myy = [0 0 2 6*y 12*y^2 20*y^3 ...
            0 0 2*x 6*x*y 12*x*y^2 ...
            0 0 2*x^2 6*x^2*y ...
            0 0 2*x^3 0 0 0]';
end

function mx = derx(x,y)
    mx = [0 0 0 0 0 0 ...
          1 y y^2 y^3 y^4 ...
          2*x 2*x*y 2*x*y^2 2*x*y^3 ...
          3*x^2 3*x^2*y 3*x^2*y^2 ...
          4*x^3 4*x^3*y 5*x^4]';
end

function my = dery(x,y)
    my = [0 1 2*y 3*y^2 4*y^3 5*y^4 ...
          0 x 2*x*y 3*x*y^2 4*x*y^3 ...
          0 x^2 2*x^2*y 3*x^2*y^2 ...
          0 x^3 2*x^3*y ...
          0 x^4 0]';
end

function z = monomials(x,y)
    z = [1 y y^2 y^3 y^4 y^5 ...
         x x*y x*y^2 x*y^3 x*y^4 ...
         x^2 x^2*y x^2*y^2 x^2*y^3 ...
         x^3 x^3*y x^3*y^2 ...
         x^4 x^4*y x^5]';
end

function M = reference()
    A(1,:)=monomials(0,0); A(2,:)=monomials(1,0); A(3,:)=monomials(0,1);
    A(4,:)=derx(0,0); A(5,:)=dery(0,0);
    A(6,:)=derx(1,0); A(7,:)=dery(1,0);
    A(8,:)=derx(0,1); A(9,:)=dery(0,1);
    A(10,:)=derxx(0,0); A(11,:)=derxy(0,0); A(12,:)=deryy(0,0);
    A(13,:)=derxx(1,0); A(14,:)=derxy(1,0); A(15,:)=deryy(1,0);
    A(16,:)=derxx(0,1); A(17,:)=derxy(0,1); A(18,:)=deryy(0,1);
    A(19,:) = dery(0.5,0);
    A(20,:) = -1*derx(0,0.5);
    A(21,:) = (-1/sqrt(2))*(derx(0.5,0.5)+dery(0.5,0.5));

I21 = eye(21);

M = zeros(21,21);
for i=1:21
    f = I21(:,i);
    M(i,:) = (A\f)';
    % If any elements are small but nonzero from matrix inversion rounding error, set to 0
    for j=1:21
        if abs(M(i,j)) < 1e-12
            M(i,j) = 0;
        end
    end
end
disp('Reference basis functions computed.');
end