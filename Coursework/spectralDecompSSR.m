function G = spectralDecompSSR(s,A,B,C,D)
% Computes the spectral decomposition of the State Space Model (A,B,C,D)
% Input: State Space Coefficient Matrices (A,B,C,D)
% Output: Resultant Transfer Function computed from the Spectral
% Decomposition/Dyadic Expansion (Eqn. 1 in Coursework)

%% Eigenvectors and Eigenvalues of the System
[right,eigVals,left] = eig(A);
% right: contains the right eigenvectors as columns
% left: contains the left eigenvectors as columns
% eigVals: diagonal matrix populated by eigenvalues 

%% Formulating the Dyadic Expansion
n = size(A,1); % # States in the system
g = NaN(n,n,n);
for ii = 1:n
    g(:,:,ii) = right(:,ii)*(left(:,ii)');
    g(:,:,ii) = g(:,:,ii)/(s - eigVals(ii,ii));
end

G = sum(g,3);
G = C*G*B + D;
end