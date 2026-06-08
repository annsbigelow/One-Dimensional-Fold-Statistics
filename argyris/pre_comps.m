%% This script defines the Argyris basis functions

A = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1;
    0 0 0 0 0 0 1 0 0 0 0 2 0 0 0 3 0 0 4 0 5;
    0 1 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 0 0 1 0;
    0 0 0 0 0 0 0 0 0 0 0 2 0 0 0 6 0 0 12 0 20;
    0 0 0 0 0 0 0 1 0 0 0 0 2 0 0 0 3 0 0 4 0;
    0 0 2 0 0 0 0 0 2 0 0 0 0 2 0 0 0 2 0 0 0;
    1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0;
    0 1 2 3 4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 2 2 2 2 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 2 3 4 0 0 0 0 0 0 0 0 0 0;
    0 0 2 6 12 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 -1 0 0 0 0 0 -1/2 0 0 0 0 -1/4 0 0 0 -1/8 0 0 -1/16 0;
    0 0 0 0 0 0 -1 -1/2 -1/4 -1/8 -1/16 0 0 0 0 0 0 0 0 0 0;
    0 1 1 3/4 1/2 5/16 1 1 3/4 1/2 5/16 1 3/4 1/2 5/16 3/4 1/2 5/16 1/2 5/16 5/16];

I21 = eye(21);

% Each column of Phi holds the coefficients defining each basis function
Phi = zeros(21,21);
for i=1:21
    f = I21(:,i);
    Phi(:,i) = A\f;
    % If any elements are small but nonzero from matrix inversion rounding error, set to 0
    for j=1:21
        if abs(Phi(j,i)) < 1e-12
            Phi(j,i) = 0;
        end
        Phi(j,i) = round(Phi(j,i),2);
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

% Careful! S and Phi use 1-based indexing
for I = 1:21
for J = 1:21
for p=0:10
for q=0:10-p
    % Define combined coefficients from basis function coefficients
    c_pq=0;
    for i=0:5
    for j=0:5-i
        m=p-i; n=q-j;
        if (m>5 || n<0 || m<0 || m+n>5) continue;
        end
        a_glbl_idx = idx(i,j);
        b_glbl_idx = idx(m,n);
        ab = Phi(a_glbl_idx,I)*Phi(b_glbl_idx,J);
        c_pq = c_pq + ab;
    
        % Force matrix
        for r=1:3
            [alpha,beta]=dmap(r);
        for s=1:3
            [delta,gamma]=dmap(s);

            cross_ab = D(alpha,2)*D(beta,1) + D(alpha,1)*D(beta,2);
            cross_dg = D(delta,2)*D(gamma,1) + D(delta,1)*D(gamma,2);
    
            a=i*m*(i-1)*(m-1)*D(alpha,1)*D(beta,1)*D(delta,1)*D(gamma,1);
            b=i*m*n*(i-1)*D(alpha,1)*D(beta,1)*cross_dg;
            c=i*n*(i-1)*(n-1)*D(alpha,1)*D(beta,1)*D(delta,2)*D(gamma,2);
            d=i*j*m*(m-1)*D(delta,1)*D(gamma,1)*cross_ab;
            e=i*j*m*n*cross_ab*cross_dg;
            f=i*j*n*(n-1)*D(delta,2)*D(gamma,2)*cross_ab;
            g=j*m*(j-1)*(m-1)*D(alpha,2)*D(beta,2)*D(delta,1)*D(gamma,1);
            h=j*m*n*(j-1)*D(alpha,2)*D(beta,2)*cross_dg;
            k=j*n*(j-1)*(n-1)*D(alpha,2)*D(beta,2)*D(delta,2)*D(gamma,2);
    
            F(I,J,r,s) = F(I,J,r,s) + ...
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
    S(I,J) = S(I,J) + c_pq*factorial(p)*factorial(q)/factorial(p+q+2);

end
end
end
end
disp('Finished assembly of stiffness and force matrices.');

%% Reshape all arrays into flat vectors for output
% Phi_out[21*i+j] = Phi_i[j], i=0,...,20, j=0,...,20
Phi_out = reshape(Phi,[1,441]);
% S_out[21*i+j] = \int_R phi_i phi_j
S_out = reshape(S,[1,441]);
% F_out[9(21*i+j) + 3r + s] = \int_R (H_i)_{ab}*(H_j)_{dg}
%   where r->(a,b) and s->(d,g) via [a,b]=dmap(r) and [d,g]=dmap(s).
% So, for each (i,j) there are 9 elements
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

%% Output all flat arrays in binary
PhiID = fopen('Phi.bin', 'w');
fwrite(PhiID, Phi_out, 'double');
fclose(PhiID);

SID = fopen('S.bin', 'w');
fwrite(SID, S_out, 'double');
fclose(SID);

FID = fopen('F.bin', 'w');
fwrite(FID, F_out, 'double');
fclose(FID);

disp('Successfully outputted flat arrays.');

%% Output flat arrays as C++ Source
cID = fopen('pre_comps.hh','w');
fprintf(cID,'extern const double S[%d];\n',21*21);
fprintf(cID,'extern const double F[%d];\n',9*21*21);
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
fprintf(cID, '};');
fclose(cID);

%% Check if the row-sums of the stiffness matrix is constant (they are not)
disp('Begin S row sums');
for i=0:20
    row_sum=0;
    for j=0:20
        row_sum = row_sum + S_out(21*i+j+1);
    end
    disp(row_sum);
end
disp('End S row sums');
for chunk=1:3
    chunk_sum=0;
    for i=1:6
        for j=1:21
            chunk_sum=chunk_sum + S(i+6*(chunk-1),j);
        end
    end
    disp(chunk_sum);
end
norm_sum=0;
for i=18:21
    for j=1:21
        norm_sum=norm_sum+S(i,j);
    end
end
disp(norm_sum);
disp('End S chunk sum');
%% Functions
% Beta function 
function z = fac(r,s)
    if r<0 || s<0
        z=0;
    else
        z = factorial(r)*factorial(s)/factorial(r+s+2);
    end
end

% Kronecker delta
function out = D(a,b)
    if a==b
        out=1;
    else
        out=0;
    end
end

% r->(alpha,beta) mapping
function [alpha,beta] = dmap(r)
    if r==1
        alpha=1; beta=1;
    elseif r==2
        alpha=1;beta=2;
    elseif r==3
        alpha=2;beta=2;
    end
end

