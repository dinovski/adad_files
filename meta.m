function a=im(people,snps, f, N, studies, p)


%people - number of people in general
%N - number of people in study
%snps - number of independent snps
%p - phenotypes


N = studies * round(N/studies);
a=(rand(people,snps)>f)+(rand(people,snps)>f); %genotypes in HW
pi = mean(a)/2; %SNP frequency in population (also called pi_star in Visscher and Hill).
phenotype = randn(people,p);

for i=1:p
    phenotype(:,p) = phenotype(:,p) - mean(phenotype(:,p));
end



for i=1:p
    %every center does a GWAS study:
    for s=1:studies
        share = N / studies;
        start = share*(s-1)+1;
        stop  = share*s;
    
        for t=1:snps
            x = a(start:stop,t);
            y = phenotype(start:stop,i);
    
            beta(t,s,i) = (x'*y)/(x'*x);
            w(t,s,i) = (sum((y-beta(t,s,i)*x).^2) / (share - 2)) / sum((x - mean(x)).^2);
            w(t,s,i) = 1/w(t,s,i);
        end
    end

    up = sum((beta(:,:,i) .* w(:,:,i)),2);
    down = sqrt(sum(w(:,:,i),2));
    Z(:,i) = up./down;
    
end






% for t=1:snps
%     beta_real(t,1) = a(1:N,t)'*phenotype(1:N)/(a(1:N,t)'*a(1:N,t));
% end



%attack 
for i=1:p
    for t=1:people
         a_ind = a(t,:);
         test(t,i) = 1/snps*((a_ind-2*pi)*Z(:,i));
    end    
end

test = abs(test);
test_big = sum(test,2);


%Draw
figure(1);clf;hold on;
plot(test_big(1:N),'or');
plot([N+1:length(test)],test_big(N+1:end),'ok');

 
figure(2);


plot(phenotype(1:N),test(1:N),'o');
figure(1);

