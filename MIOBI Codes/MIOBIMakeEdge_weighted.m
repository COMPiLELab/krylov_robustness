function [robScoreNoUpdate, robScoreUpdate50, Aorg] = MIOBIMakeEdge_weighted(fname,E,k,topT)
% filename: in format: source dest (integer)
% E: search space of weighted edges
% k: budget
% top t: top t eigenvalue and eigenvector to compute the rob. score
% robScoreUpdate50 - recompute at every 50 operations and update using
% perturbation theory
% robScoreNoUpdate - no recompute at all 

% set the number of eigenpairs use in the perturbatin update 
t = 50;

format longG

[a,b,~] = find(fname);

data = [a b];

N = max(max(data(:,1)),max(data(:,2)));

%K = [10 50 100 200 350];
%K = [0.01*N 0.05*N 0.1*N 0.15*N 0.25*N];
%K = round(K)
%k = max(K);

                    
% build the matrix
A = sparse(data(:,1), data(:,2), data(:,3), N, N);
Aorg = A;

if ~ishermitian(A)
	error('MIOBIBreakEdge2_weighted:: adjacency matrix has to be symmetric');
end

        
R = zeros(1+1,1);
Rtrue = zeros(1+1,1);

en = eigs(A,topT);
        
R(1) = log((1/topT)*sum(exp(en)));
Rtrue(1) = log((1/topT)*sum(exp(en)));




[V,D] = eigs(A, t);
V = [abs(V(:,1)) V(:,2:t)];


% returning the eignvalue
e = diag(D); %D(D~=0);



for j = 1:k
    
    %degree = sum(Aorg,2);
    %maxDegree = max(degree);
    %find k+dmax with hihgest u1 eigen-scores
    %[~,sortIndex] = sort(V(:,1),'descend');
    %maxIndex = sortIndex(1:(min(maxDegree+k,N)));
    
    %pAdd = zeros((maxDegree+k),1);
    %rAdd = zeros((maxDegree+k),1);
    
    %numE = 0;
    numE = size(E, 1);
    %for n=1:length(maxIndex)
    %    for m=n+1:length(maxIndex)
    %        if Aorg(maxIndex(n),maxIndex(m)) == 0
    %            numE = numE + 1;
    %            pAdd(numE) = maxIndex(n);
    %            rAdd(numE) = maxIndex(m);
    %        end
    %    end
    %end
	pAdd = E(:, 1);
    rAdd = E(:, 2);
    pAdd(pAdd==0) = [];
    rAdd(rAdd==0) = [];
    
    %procedure starts here
    Up = V(pAdd,:); %getting the eigenvector value at p
    Ur = V(rAdd,:); %getting the eignevector value at r
    c = exp(e)'; %-e(1))*exp(e(1)); %computing the constant c
    
    c = c(ones(1,length(pAdd)),:);
    
    compScore = times(c,exp(2.*times(Up,Ur))); %computing score for each (sum) component in Equ (8)
    
    score = sum(compScore,2);
    
    %returning the index of all the lowest score edges
    [maxscore,maxIndex] = max(score);
    %minIndex = find(score==min(score));
    
    %pick an arbitrary one
    add = [pAdd(maxIndex(1)) rAdd(maxIndex(1))];
    
    %S
    %pause
    
    Aorg(add(1), add(2)) = E(maxIndex(1), 3);
    Aorg(add(2), add(1)) = E(maxIndex(1), 3);
    
    deltaA = sparse(N,N);
    %update top t eigenvalues of A equation (4)
    deltaA(add(1), add(2)) = E(maxIndex(1), 3);
    deltaA(add(2), add(1)) = E(maxIndex(1), 3);

	E = E([1:maxIndex(1)-1, maxIndex(1)+1:numE], :); % remove edge from the search space
    
    deltaH = (V'*deltaA*V); %
    
    deltaE = diag(deltaH); %diagonal is the values we want
    
    tildaE = e + deltaE;
    
    deltaH(logical(eye(size(deltaH)))) = 0;
    
    %(V(:,2)'*deltaA*V(:,1))/2
    %(V(:,3)'*deltaA*V(:,1))/2
    %(V(:,4)'*deltaA*V(:,1))/2
    %(V(:,5)'*deltaA*V(:,1))/2
    %(V(:,6)'*deltaA*V(:,1))/2
    
    %update top t eigenvectors of A equation (9)
    diffE = e';
    diffE = diffE(ones(1,t),:) - e(:,ones(1,t)); % compute lambda_j - lambda_i
    diffE = rdivide(1,diffE); % compute 1/(lambda_j-lambda_i)
    %diffE(logical(eye(size(diffE)))) = 0; %set all the diagonal to zero
    [naR, naC] = find(isfinite(diffE));
    diffE(naR, naC) = 0;
    
    
    
    %diffE = abs(diffE);
    
    innerSum = times(deltaH, diffE); %the corresponding entry multipication gives the constant terms inside the sum
    
    % second way of updating the eigenvectors using loops
    %    deltaV = [];
    
    %     for i=1:t
    %         deltaVTemp = 1/2 * V(:,i);
    %         %pair = [];
    %         for k=1:t
    %             if i ~= k
    %                 %pair = [pair; i k (V(:,k)'*deltaA*V(:,i))  e(i)-e(k) (V(:,k)'*deltaA*V(:,i))/(e(i)-e(k))];
    %                 %((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k)
    %                 deltaVTemp = deltaVTemp + ((V(:,k)'*deltaA*V(:,i))/((e(i)-e(k))))*V(:,k);
    %             end
    %         end
    %         %pair
    %         deltaV = [deltaV deltaVTemp];
    %     end
    
    %    V + deltaV
    
    % another way of updating the eigenvectors
    deltaV = zeros(N,t);
    
    for i=1:t
        deltaVTemp = innerSum(:,i)';
        deltaVTemp = deltaVTemp(ones(length(V),1),:);
        %times(deltaVTemp, V)
        temp = sum(times(deltaVTemp, V),2);
        deltaV(:,i) = temp;
    end
    %V - deltaV
    
    tildaV = V + deltaV;
    
    for i=1:size(tildaV,2)
        tildaV(:,i) = tildaV(:,i)/norm(tildaV(:,i));
    end
    
    e = tildaE;
    V = tildaV;
    
    V = [abs(V(:,1)), V(:,2:t)];
    
end

enTrue = eigs(Aorg,topT);
Rtrue(2) = log((1/topT)*sum(exp(enTrue)));
robScoreNoUpdate = ((Rtrue(2)-Rtrue(1))*100/Rtrue(1));

end

