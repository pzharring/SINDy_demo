function Xi = SparseRegression(Theta, dXdt, lambda, n)
% Perform the sparse regression 
% INPUTS:  Theta  -- matrix of potential RHS functions
%          dXdt   -- matrix of the X dot data
%          lambda -- sparsification parameter
%          n      -- dimensionality of the problem (number of vars in x)
% OUTPUTS: Xi     -- Vector of coefficients for the RHS functions

Xi = Theta\dXdt;
for k=1:10
   smallinds = (abs(Xi)<lambda);
   Xi(smallinds)=0;
   for ind = 1:n 
      biginds = ~smallinds(:,ind);
      Xi(biginds,ind) = Theta(:,biginds)\dXdt(:,ind);
   end
end