function a=im(people,snps, f, N)

%people - number of people in general
%N - number of people in study
%snps - number of independent snps

a=(rand(people,snps)>f)+(rand(people,snps)>f); %genotypes in HW
pi = mean(a)/2; %SNP frequency in population (also called pi_star in Visscher and Hill).
p_sample = mean(a(1:N,:))/2;

phenotype = randn(people,1);

for t=1:snps
    beta(t,1) = a(1:N,t)'*phenotype(1:N)/(a(1:N,t)'*a(1:N,t));
    
end

beta_sign = sign(beta);

for t=1:people
    a_ind = a(t,:);
    test(t) = N/snps * ((a_ind-2*pi)*beta_sign);
end

%Draw
figure(1);clf;hold on;
plot(test(1:N),'or');
plot([N+1:length(test)],test(N+1:end),'ok');


figure(2);


plot(phenotype(1:N),test(1:N),'o');
    



x=1;

