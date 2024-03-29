function [allA,allB,allC,allVarExps]=bootstrappedPARAFAC(X, individual_metadata, numFactors, maxIterations)
[I,J,K] = size(X);

Options = [];
Options(1) = 1e-6;  % convergence (default 1e-6)
Options(2) = 2;     % initialization (default 1)
Options(3) = 0;     % resulting plot (default 0)
Options(4) = 2;     % no scaling (default 1)
Options(5) = NaN;	% no output (default 10)
Options(6) = 2500;   % max number of iterations (default 2500)

const = [0 0 0];    % PARAFAC mode constraints, currently none

% Find individual groups
RFlow = find(individual_metadata(:,2) == "0");
RFmid = find(individual_metadata(:,2) == "1");
RFhigh = find(individual_metadata(:,2) == "2");

% Establish which individuals to remove randomly
combinations = nchoosek(1:I, 3);
combinations_low = ismember(combinations, RFlow);
combinations_mid = ismember(combinations, RFmid);
combinations_high = ismember(combinations, RFhigh);
valid = sum(combinations_low,2) == 1 & sum(combinations_mid,2) == 1 & sum(combinations_high,2) == 1;
combinations = combinations(valid,:);
%maxIterations = size(combinations, 1);

if (maxIterations == size(combinations, 1))
    jackKnifedSamples = 1:maxIterations;
else
    jackKnifedSamples = randi([1 size(combinations,1)],1,maxIterations);
end

allA = [];
allB = [];
allC = [];
allVarExps = zeros(size(jackKnifedSamples,2),1);
textprogressbar('Creating bootstrapped PARAFAC models ');

for n=1:size(jackKnifedSamples,2)
    textprogressbar((n/size(jackKnifedSamples,2))*100);
    Xboot = X;
    %IDremoval = [RFlow(RFlow_remove(n)) RFmid(RFmid_remove(n)) RFhigh(RFhigh_remove(n))];
    selection = jackKnifedSamples(n);
    IDremoval = combinations(selection,:);
    Xboot(IDremoval,:,:) = [];
    
    evalc('[Factors, it, err, c] = silent_parafac(Xboot, numFactors, Options, const)');

    [A, B, C] = fac2let(Factors);
    [A, B, C] = sortParafacComponents(Xboot, A, B, C);

    % Make ID loading vector with NaNs for the removed individuals
    Amodified = zeros(I,numFactors);
    IDiterator = 1;
    for i=1:I
        if(ismember(i,IDremoval))
            Amodified(i,:) = nan(1,numFactors);
        else
            Amodified(i,:) = A(IDiterator,:);
            IDiterator = IDiterator + 1;
        end
    end

    allA(:,:,n) = Amodified;
    allB(:,:,n) = B;
    allC(:,:,n) = C;
    allVarExps(n) = calcVarExplained(Xboot, A, B, C);
end

evalc('[Factors, it, err, c] = silent_parafac(X, numFactors, Options, const)');
[realA, realB, realC] = fac2let(Factors);
[allA, allB, allC] = checkComponentNegativity(realA, realB, realC, allA, allB, allC);

textprogressbar(' done');