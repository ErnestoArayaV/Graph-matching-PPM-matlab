%Algorithm for the initialization of the PPM method at random on a
%Frobenius neighborhood of the ground truth permutation. Use a derrangement
%generator randpermfull below.
%Input:
%n               <----- size of the permutation
%initial_ball    <----- number of fixed points of the output permutation
%P_0             <----- grounth truth permutation
%Outout:
%P_init          <----- a random permutation matrix within the specified Frobenius
%distance of P_0
function P_init=initial_perm(n,initial_ball,P_0)
    P_init=P_0;                             %latent (g.t) permutation
   %P_rnd = eye(n);                        %if we assume g.t=Id
    p1=randperm(n,initial_ball);           %choose 'initial_ball' indices at random
    p2=p1(randpermfull(length(p1)));       %generate a derrangement on those indices
    P_init(:,p2) = P_init(:,p1);