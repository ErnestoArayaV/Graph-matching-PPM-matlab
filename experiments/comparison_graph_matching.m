%PPM for seeded graph matching 
%Numerical experiments for the article: https://arxiv.org/abs/2204.04099 


%% initialization
%p = 0.5;%only for the methods with community structure
%q=0.2;
%d=3;% only for the methods with geometry...Wishart, RGG, etc...
frac_of_fps=0.01; %fraction of initial fixed points(seed accuracy, for PPM initialization....
maxiter = 7; %maximum number of projected power iterations 
num_run = 25;%number of Montecarlo runs 
vec_dim =  800;%500:500:1500;%50:200:850;%size of the matrix (can be an array)
len_dim = length(vec_dim);
vec_noise =0.5:0.05:1;%'noise' parameter sigma (can be an array)
len_noise = length(vec_noise);
pi_corr = zeros(len_dim, len_noise, num_run);% performance measure of the estimator1
pi2_corr = zeros(len_dim, len_noise, num_run);% performance measure of the estimator2
pi3_corr = zeros(len_dim, len_noise, num_run);% idem
pi4_corr = zeros(len_dim, len_noise, num_run);%
tic
%% Montecarlo runs
for ind_run = 1:num_run
    %% Iterate over dimensions
    for ind_dim = 1:len_dim
        n = vec_dim(ind_dim);
        %fraction of fps and initialization 
        in_ball=floor((1-frac_of_fps)*n); 
        %% Iterate over noise levels 
        for ind_noise = 1:len_noise
            sigma = vec_noise(ind_noise);
            %Generate random A,B according to a model
            [A, B, P_rnd] = generate_wig(n,sigma);%A,B from the correlated Wigner model
            %[A, B, A0, B0, P_rnd] = generate_er(n, p, sigma);%A,B from the correlated Erdos-Renyi model
            %c1=n/2;c2=n/2;
            %[A,B,A0,B0,P_rnd]=generate_SMB_2(c1,c2, p,q,sigma);%A,B from the correlated SBM model
            %[A, B, P_rnd] = generate_wish(n,d,1,sigma);%A,B from the correlated Wishart model
           
            %% Projected power method      
            %define the initial permutation from the seed
            P_init =initial_perm(n,in_ball,P_rnd); %%%%%%Initialization matrix
            P = matching_proj_it_in(A, B,maxiter,P_init);
            %% Grampa method 
            P2= matching_robust_spectral(A,B,0.2);
            %% Umeyama method
            P3= matching_umeyama(A,B);
            %% Top eigenvector method
            P4= matching_top_eigvec(A,B);
            %% Compute performance measures
            fix_pt_ratio = sum(dot(P_rnd, P)) / n;
            fix_pt_ratio2= sum(dot(P_rnd, P2)) / n;
            fix_pt_ratio3= sum(dot(P_rnd, P3)) / n;
            fix_pt_ratio4= sum(dot(P_rnd, P4)) / n;
            pi_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio;
            pi2_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio2;
            pi3_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio3;
            pi4_corr(ind_dim, ind_noise, ind_run) = fix_pt_ratio4;
            toc;
        end
    end
end
%comput the mean for the performance measure over all montecarlo runs
m_pi_corr=mean(pi_corr,3);
m_pi2_corr=mean(pi2_corr,3);
m_pi3_corr=mean(pi3_corr,3);
m_pi4_corr=mean(pi4_corr,3);
%generate plot for the performance measure
for j=1:len_dim
    figure;hold on;
    hdata1=line(vec_noise, m_pi_corr(j,:));
    hdata2=line(vec_noise, m_pi2_corr(j,:));
    hdata3=line(vec_noise, m_pi3_corr(j,:));
    hdata4=line(vec_noise, m_pi4_corr(j,:));
    set(hdata1, 'LineStyle','--','Marker', 'o', 'MarkerSize', 5,'MarkerEdgeColor', 'none', 'MarkerFaceColor', [.75 .75 1],'LineWidth',2.0);
    set(hdata2, 'LineStyle','--','Color', 'r','Marker', 'x', 'MarkerSize', 5,'MarkerEdgeColor', 'r', 'MarkerFaceColor', [.75 .75 1],'LineWidth',2.0);
    set(hdata3, 'LineStyle','--','Color', [0 0.5 0],'Marker', '*', 'MarkerSize', 5,'MarkerEdgeColor', [0 0.5 0], 'MarkerFaceColor', [.75 .75 1],'LineWidth',2.0);
    set(hdata4, 'LineStyle','--','Color', [.5 0 .5],'Marker', '^', 'MarkerSize', 5,'MarkerEdgeColor', [.5 0 .5], 'MarkerFaceColor', [.75 .75 1],'LineWidth',2.0);
    str = sprintf('Compare matching methods for n=%i, Wigner model',vec_dim(j));
    hTitle=title(str);set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')
    set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02],'XMinorTick', 'on', 'YMinorTick', 'on','XColor', [.4 .4 .4], 'YColor', [.4 .4 .4],'LineWidth', 1,'Fontsize',14);
    %plot(vec_noise, m_pi_corr(j,:));hold on;plot(vec_noise, m_pi2_corr(j,:));hold on;plot(vec_noise, m_pi3_corr(j,:));hold on;plot(vec_noise, m_pi4_corr(j,:));
    hXLabel = xlabel('Noise $\sigma$','interpreter','latex');
    hYLabel = ylabel('Recovery fraction');
    set([hXLabel, hYLabel], 'FontSize', 20);
    hLegend = legend([hdata1, hdata2,hdata3,hdata4], 'PPM in.1','PPM in.2','PPM in.3','PPM in.4');%,'Top-eigen','Top-eigen+PPM');
    set(hLegend, 'FontSize', 11);hold off
    
end

           
