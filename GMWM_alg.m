%Greedy algorithm GMWM, maximum version
%algorithm introduced in "Joseph Lubars and R. Srikant. Correcting the output of approximate graph matching algorithms. IEEE
%conference on Computer Communications, pages 1745–1753, 2018."...
%....this version is the one with erase step described in "Seeded graph
%matching for the correlated Wigner model via the projected power method"
%arXiv:2204.04099
%Inputs: 
%M           <------ matrix of 'costs' 
%param       <------ negative parameter to 'emulate' the erase step
%Output:
%X           <------ a permutation matrix 

function X=GMWM_alg(M,param)
    n=length(M);
    X=zeros(n,n);
    for i=1:n
        [~,I] = max(M(:));
        [I_row, I_col] = ind2sub(size(M),I);
        X(I_row, I_col)=1;
        M(I_row,:)=param*ones(1,n);
        M(:,I_col)=param*ones(n,1);
    end 
        