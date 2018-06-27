function [Choice, Offers] = choicemat(S, varargin)
%this provides the choice matrix with constraints entered. S is the
%structure containing the data. constraints should be entered as strings
%followed by a boolean of how to include the constraint 
%i.e. if laser only trials are used, put in ('lasert', true)
%Choice_MAT is output as A choices on the first row and B choices on the
%second row
use = true(size(S.offer));
Offers = S.Offers;
for i = 1:2:numel(varargin)
    use = use & noty(varargin{i + 1}, S.(varargin{i}));
end    
%Choice = fliplr([sum(S.offtype(use & S.offer,:)); sum(S.offtype(use & ~S.offer,:))]);

Choice =([sum(S.offtype(use & S.offer,:)); sum(S.offtype(use & ~S.offer,:))]);
end

function [X] = noty(y,x)
if y
    X = x; 
else
   X = ~x;
end 
end