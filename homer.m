
function a=homer(people,snps, f, N,top)

%people - number of people in general
%N - number of people in study
%snps - number of independent snps


a=(rand(people,snps)>f)+(rand(people,snps)>f); %genotypes in HW
pi = repmat((1-f),1,snps);%mean(a)/2; %SNP frequency in population (also called pi_star in Visscher and Hill).
p_sample = mean(a(1:N,:))/2;


%GWAS
% cases=(rand(N,snps)>f)+(rand(N,snps)>f);
% p_case = mean(cases(1:N,:))/2;
% or = (p_sample./p_case) ./ ((1-p_sample)./(1-p_case));
% or=abs(log(or));
% %or=rand(1,snps);
% [ b, ix ] = sort( or(:), 'descend' );
% ix_top = ix(1:top);
% p_sample = p_sample(ix_top);
% pi = pi(ix_top);
% snps = top;
% a = a(:,ix_top);
% %



vec_left  = log(1 - p_sample+eps) - log(1 - pi);
vec_right = log(p_sample+eps) - log(pi);

for t=1:people
    a_ind = a(t,:);
    g0 = (a_ind == 0);
    g1 = (a_ind == 1);
    g2 = (a_ind == 2);
    
    test(t) = (2*g0 + g1) * vec_left' + (g1 + 2*g2) * vec_right';
end

expected_value = 0.5*snps/N;

%Draw
figure(1);clf;hold on;
plot(test(1:N),'or');
plot([N+1:length(test)],test(N+1:end),'ok');

line([1 N/2],[expected_value expected_value],'Color',[0.5 0 0],'LineWidth',3); %execptation for positive tests based theoretical estimation
line([N/2 N],[mean(test(1:N)) mean(test(1:N))],'Color',[0.5 0 0],'LineStyle','--','LineWidth',3) %what we saw in the test group for positives

line([N+1 N+1+(people-N)/2],[-expected_value -expected_value],'Color',[0.5 0.5 0.5],'LineWidth',3) %expectation for negative tests based on theoretical estimation
line([N+1+(people-N)/2 people],[mean(test(N+1:people)) mean(test(N+1:people))],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',3) %expectation for negative tests based on theoretical estimation





x=1;


test_sorted = sort(test);
true_test = test(1:N);
false_test = test(N+1:length(test));

for t=1:length(test_sorted)
    thresh = test_sorted(t);
    sensitivity(t) = length(find(true_test >= thresh)) / length(true_test);
    specificity(t) = length(find(false_test < thresh)) / length(false_test);
end

figure(2);
plot(log10(1-specificity),sensitivity,'o-');

    
    



