function a = wang(left, right, common_freq, study, test)

%get rid of rare SNP:
common_freq = 1-common_freq;
common = find(0.5*(mean(left')+mean(right'))<common_freq);
left_common = left(common,:);
right_common = right(common,:);


%divde the group to a study and a reference panel:
pos = randperm(size(left,2));

cases = pos(1:study); %people who participated in the study
negative = pos(study+1:study+1+test); %people who did not participate
control = pos(study+1+test+1:end); %people who did not participate

left_cases = left_common(:,cases);
right_cases = right_common(:,cases);

left_negative = left_common(:,negative);
right_negative = right_common(:,negative);

left_control = left_common(:,control);
right_control = right_common(:,control);





%estimate 
r_reference = 0.5*(corr(left_control') + corr(right_control'));
r_cases = 0.5*(corr(left_cases') + corr(right_cases'));

delta = r_cases - r_reference;
for t=1:size(delta,2)
    delta(t,t) = 0;% a bit ugly.
end

imagesc(delta);

% figure(1);
% imagesc(r_reference.^2);
% figure(2);
% imagesc(r_cases.^2);

for t=1:length(cases)
    victim_left = 2*left_cases(:,t)-1;
    victim_right = 2*right_cases(:,t)-1;
    victim_corr = victim_left * victim_left' + victim_right * victim_right';
    Tr_cases(t) = sum(sum(delta.*victim_corr))/2;
end

for t=1:length(negative)
    victim_left = 2*left_negative(:,t)-1;
    victim_right = 2*right_negative(:,t)-1;
    victim_corr = victim_left * victim_left' + victim_right * victim_right';
    Tr_negative(t) = sum(sum(delta.*victim_corr))/2;
end



clf;
plot(Tr_cases,'o-r')
hold on
plot([study+1:study+1+test], Tr_negative,'o-k')
a=Tr;