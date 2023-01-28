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
m = size(B,2); % # Inputs in the system
l = size(C,1); % # Outputs in the system
u = NaN(m,1,n);
y = NaN(l,1,n);
G = D;
for ii = 1:n
    g = NaN;
    u(:,:,ii) = B'*left(:,ii); % Input Pole Vector
    y(:,:,ii) = C*right(:,ii); % Output Pole Vector
    lMag(ii) = norm(left(:,ii),2); % Magnitude (2-Norm) of left eigenvectors
    rMag(ii) = norm(right(:,ii),2); % Magnitude (2-Norm) of right eigenvectors
    dot(ii) = left(:,ii)'*right(:,ii);
    g = (y(:,:,ii)*u(:,:,ii)')/(s - eigVals(ii,ii));
    G = G + g;
end

end