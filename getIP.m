function [IP, Sigma] = getIP(ChoiceMAT, Offers)

x = log(Offers(1, 2:end- 1)./Offers(2, 2:end- 1));
%figure
%scatter(1:numel(x), Choice(:,2)./sum(Choice, 2), 'k')
Choice(1,:) = ChoiceMAT(2,2:end - 1);
Choice(2,:) = sum(ChoiceMAT(:,2:end - 1));

[b] = glmfit(x, Choice', 'binomial', 'link', 'probit');
IP = -b(1)/b(2);
Sigma = -1/b(2);
