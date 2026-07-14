%% This script defines the Argyris basis functions
clear;
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

% Each row of M holds the coefficients defining each basis function
M = zeros(21,21);
for i=1:21
    f = I21(:,i);
    M(i,:) = (A\f)';
    % If any elements are small but nonzero from matrix inversion rounding error, set to 0
    for j=1:21
        if abs(M(j,i)) < 1e-12
            M(j,i) = 0;
        end
    end
end
disp('Basis functions computed.');
%% Now that we have the coeffs for each basis function, fill the local stiffness matrix,
% whose entries are the integrated pairwise products of all the basis
% functions.

% Also compute the integrated Laplacian products.
% For each pairwise product (there are 21x21 of these), there will be a
% rank-4 tensor, depending on the variables that are differentiated with
% respect to.

S = zeros(21,21); % Local stiffness matrix
% Local forces: huge tensor
% F is indexed as (I,J,alpha,beta,delta,gamma)
F = zeros(21,21,3,3); 
idx = @(i,j) 1 + 6*i - i*(i-1)/2 + j;

% Careful! S and M use 1-based indexing
for I = 1:21
for J = 1:21
    for i=0:5
    for j=0:5-i
    for m=0:5
    for n=0:5-m
        p=i+m; q=j+n;
        a_glbl_idx = idx(i,j);
        b_glbl_idx = idx(m,n);
        ab = M(I,a_glbl_idx)*M(J,b_glbl_idx);

        S(I,J) = S(I,J) + ab*fac(p,q);
    end
    end
    end
    end
end
end
disp(cond(S));
writematrix(S,"S_precomp.csv");
%%
% Force tensor
for I=1:21
for J=1:21
    for i=0:5
    for j=0:5-i
    for m=0:5
    for n=0:5-m
        a_glbl_idx = idx(i,j);
        b_glbl_idx = idx(m,n);
        ab = M(I,a_glbl_idx)*M(J,b_glbl_idx);
        p=i+m; q=j+n;
        for r=0:2
        for s=0:2
            a=i*m*(i-1)*(m-1)*delta(r,0)*delta(s,0);
            b=i*m*n*(i-1)*delta(r,0)*delta(s,1);
            c=i*n*(i-1)*(n-1)*delta(r,0)*delta(s,2);
            d=i*j*m*(m-1)*delta(r,1)*delta(s,0);
            e=i*j*m*n*delta(r,1)*delta(s,1);
            f=i*j*n*(n-1)*delta(r,1)*delta(s,2);
            g=j*m*(j-1)*(m-1)*delta(r,2)*delta(s,0);
            h=j*m*n*(j-1)*delta(r,2)*delta(s,1);
            k=j*n*(j-1)*(n-1)*delta(r,2)*delta(s,2);
    
            F(I,J,r+1,s+1) = F(I,J,r+1,s+1) + ...
                ab*( ...
                a*fac(p-4,q) ...
                +(b+d)*fac(p-3,q-1) ...
                +(c+e+g)*fac(p-2,q-2) ...
                +(f+h)*fac(p-1,q-3) ...
                +k*fac(p,q-4) ...
                );
        end
        end
    end
    end
    end
    end
end
end

disp('Finished assembly of stiffness and force matrices.');

%% Reshape all arrays into flat vectors for output
% Phi_out[21*i+j] = Phi_i[j], i=0,...,20, j=0,...,20
%Phi_out = reshape(Phi,[1,441]);
M_out = reshape(permute(M,[2,1]),[],1);
% S_out[21*i+j] = \int_R phi_i phi_j
S_out = reshape(S,[1,441]);
% F_out[9(21*i+j) + 3r + s] = \int_R (H_i)_{ab}*(H_j)_{dg}
%   where r->(a,b) and s->(d,g) via [a,b]=dmap(r) and [d,g]=dmap(s).
% So for each (i,j) there are 9 elements.
% r,s=0,1,2.
F_out = reshape(permute(F,[4 3 2 1]), [], 1);

% Check if F was flattened correctly
for I=0:20
for J=0:20
for r=0:2
    for s=0:2
        if F(I+1,J+1,r+1,s+1) ~= F_out(9*(21*I+J)+3*r+s+1)
            disp('F was flattened incorrectly.');
        end
    end
end
end
end

disp('Finished flattening arrays.');

%% Output flat arrays as C++ Source
cID = fopen('pre_comps.hh','w');
fprintf(cID,'extern const double S[%d];\n',21*21);
fprintf(cID,'extern const double F[%d];\n',9*21*21);
fprintf(cID,'extern const double M[%d];\n',21*21);
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

fprintf(cID, 'const double F[%d] = {\n',9*21*21);
for k = 1:9*21*21
    fprintf(cID, '%.17g', F_out(k));
    if k < 9*21*21
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
% Beta function 
function z = fac(r,s)
    if r<0 || s<0
        z=0;
    else
        z = factorial(r)*factorial(s)/factorial(r+s+2);
    end
end

% Kronecker delta
function out = delta(a,b)
    if a==b
        out=1;
    else
        out=0;
    end
end

