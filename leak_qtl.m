function [r2, a, y, m, y_hat] = leak_qtl(people, snps, f)

a=(rand(people,snps)>f)+(rand(people,snps)>f); %genotypes
y=randn(people,1); %phenotypes

%GWAS
for t=1:snps
    m(t,1) = a(:,t)'*y/(a(:,t)'*a(:,t));
end

%Leak:
avg_phenotype_hat = mean(m)*2*(1-f)*snps;
y_hat = a * m - avg_phenotype_hat;

plot(y,y_hat,'o');
r2 = regstats(y,y_hat,'linear','rsquare');



