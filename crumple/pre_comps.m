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
    % If any elements are small, set to 0
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
F = zeros(21,21,2,2,2,2); 
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
        for alpha=1:2
        for beta=1:2
        for delta=1:2
        for gamma=1:2
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
    
            F(I,J,alpha,beta,delta,gamma) = F(I,J,alpha,beta,delta,gamma) + ...
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
% F_out[16(21*i+j) + 8a + 4b + 2d + g] = 
%                               \int_R (H_i)_{ab}*(H_j)_{dg}
% So, for each (i,j) there are 16 elements, where gamma (g) is the 
%   fastest index, delta (d) is the next fastest, then beta (b), then 
%   alpha (a).
% a,b,d,g=0,1.
F_out = reshape(permute(F,[6 5 4 3 2 1]), [], 1);
disp('Finished flattening arrays.');
%% Output all flat arrays
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

%%
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


