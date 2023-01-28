function G = jordonDecompSSR(A)
% Computes the spectral decomposition of the State Space Model (A,B,C,D)
% Input: State Space Coefficient Matrices (A,B,C,D)
% Output: Resultant Transfer Function computed from the Spectral
% Decomposition/Dyadic Expansion (Eqn. 1 in Coursework)

%% Eigenvectors and Eigenvalues of the System
[right,eigVals,left] = eig(A);
% right: contains the right eigenvectors as columns
% left: contains the left eigenvectors as columns
% eigVals: diagonal matrix populated by eigenvalues 

%% Eigenvalue Multiplicity
[uniqueEigs,~,ic] = unique(eig(A)); % Unique Eigenvalues
count = accumarray(ic,1);
mul = [uniqueEigs,count];

end