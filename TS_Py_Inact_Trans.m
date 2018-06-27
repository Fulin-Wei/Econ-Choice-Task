function [Inact,Trans] = TS_Py_Inact_Trans(TS)
%This function determines the differences in behavior between the Laser on and no laser trials
%during inactivation sessions

%For sessions measuring trnsititivity, the sessions are named as followed
%in the google spreadsheet and in the TS struct:

%Hal50 P# D# T# - # represents a number

%The P refers to the current pair, the D refers to the number of
%consecutive days the animal has received that particular pair and T refers
%to the transitivity measurement. A rat must receive 3 pairs for 1
%transitivity point
%This sets the string number to look for each of the variables
P = 11;
D = 8;
T = 14;


Inact = struct;
Trans = struct;
%This loop runs through the TS struct to create a data matrix.
nt = 0;
ni = 0;
for i = 1:numel(TS)
    
    if isempty(TS(i).Date)
        continue
    end
    
    %Inact50_Exps = strncmp('HAL50', TS(i).Experiment, 5);
    %trans = sort(unique(cellfun(@(x) str2double(x(T)), TS(i).Experiment(Inact50_Exps))));
    
    %THis determines the number of transitivity sets which were run for
    %each animal
    trans = sort(unique(cellfun(@(x) str2double(x(T)), TS(i).Experiment)));

    %t keeps track of the transitivity number
    for t = trans
        IP = [];
        offer = cell(1,3);
        trans_sessions = double(3);
        %This goes through each pair
        for p = 1:3
            
            %Finally this steps through each day
            for d = 1:3
                
                %This finds the session fitting the t, p, and d
                %curr_sess = cellfun(@(x) all([str2double(x(D)) str2double(x(P)) str2double(x(T))]  == [d p t]), TS(i).Experiment(Inact50_Exps));
                curr_sess = cellfun(@(x) all([str2double(x(D)) str2double(x(P)) str2double(x(T))]  == [d p t]), TS(i).Experiment);
                
                %If there are multiple sessions that fit the description, use
                %the other session if this is a NaN. There shouldn't be more
                %than two, hence no while loop
                sess_ind = find(curr_sess);
                
                while numel(sess_ind) > 1
                    
                    display(strcat('!!!!!There are multiple sessions listed for rat ', TS(i).Rat, ' Day ', num2str(d), ', Pair ', num2str(p) , ', Run ', num2str(t)));
  
                    if sum(TS(i).Data(curr_sess(1)).offertype(:)) < 66
                        
                        sess_ind(1) = [];
                        
                    else     
                        sess_ind = sess_ind(1);
                    end
                end
                
                if ~any(sess_ind)
                    
                    display(strcat('!!!!!There is not session listed for rat ', num2str(TS(i).Rat), ' Day ', num2str(d), ', Pair ', num2str(p) , ', Run ', num2str(t)));
                    
                    trans_sessions(p,d) = NaN;
                    Inact.IP(ni,1:2,:) = NaN;
                    Inact.Sigma(ni,1:2,:) = NaN;
                    Inact.RT(ni,1:2,:) = NaN;
                    Inact.Omissions(ni,1:2,:) = NaN;
                    Inact.Aborted(ni,1:2,:) = NaN;
                    Inact.Trials(ni,1:2,:) = NaN;
                    continue
                    
                end    

                trans_sessions(p,d) = sess_ind;
                  
            end
            
            ni = ni + 1;
            [Inact.IP(ni,:,:), Inact.Sigma(ni,:,:), Inact.RT(ni,:,:),Inact.ID(ni,:), Inact.Omissions(ni,:,:), Inact.Aborted(ni,:,:), Inact.Trials(ni,:,:)] = Inact_Data(TS(i),trans_sessions(p,:));
            
            
        end
        nt = nt + 1;
        
        if any(isnan(trans_sessions(:)))
            IP(nt,:,:,:) = NaN;
            Sigma(nt,:,:,:) = NaN;
        else
   
        [Trans.IP(nt,:,:,:), Trans.Sigma(nt,:,:,:)] = Trans_Data(TS(i),trans_sessions);
        end
        
    end                   
end

end                  

function [IP, Sigma, RT,Identitifier,Omissions,Aborted,Trials] = Inact_Data(TS,inact_sessions)
%This function determines the IP for inactivation sessions
%All IPs need to be flipped in one direction for the analysis.
%To do this, we clump all sham, opto, laser on and laser off
%trials together, then compute the IP to determine the average
%preference across all four trial types

%First determine which offers were used in each session. If the
%offers were different, combine the sessions
for n = 2:3
    
    Sess_Offers{n - 1} = TS.Data(inact_sessions(n)).Offers;
    
end

if any(Sess_Offers{1}(:) ~= Sess_Offers{2}(:))
    %THis occurs only if the offers don't match up from day to day.
    disp('Two sessions dont have the same offers.')
    IP = nan(2);
    Sigma = nan(2);
    RT = nan(2);
    Omissions = nan(2);
    Aborted = nan(2);
    Trials = nan(2);
    Identitifier = {strcat(num2str(TS.Rat),'-',TS.Date{inact_sessions(n)}),[]};
    return
else
end

All_Choice = zeros(size(Sess_Offers{1}));
for n = 2:3

    All_Choice = All_Choice + choicemat(TS.Data(inact_sessions(n)));
    CableType{n} = TS.Inact{inact_sessions(n)};
end

[IPAll, ~] = ModelFit(All_Choice, TS.Data(inact_sessions(n)).Offers);

%flip is a boolean determining whether to flip the IP. Note that Sigma does
%not get flipped
flip = IPAll < 0;

%Inact.IP and Inact.Sigma are 3D matrices with a 2 day session of sham and opto as the first dimension,
%sham/opto as the second dimension, and laser off/on as the
%third dimension

%This checks whether this was a sham or opto session
for n = 2:3
    
    if any(strcmpi(CableType{n}, 'SHAM'))
    
        [inacttemp,Sigma(1,:)] = ModelFitStruct(TS.Data(inact_sessions(n)));
        if flip
            IP(1,:) = -1*inacttemp;
        else
            IP(1,:) = inacttemp;
        end
        %This gets only reaction times for non-nan trials and for non-laser
        %trials
        RT(1,1) = mean(TS.Data(inact_sessions(n)).rts(~TS.Data(inact_sessions(n)).lasert & ~isnan(TS.Data(inact_sessions(n)).rts)));
        
        %This is for laser
        RT(1,2) = mean(TS.Data(inact_sessions(n)).rts(TS.Data(inact_sessions(n)).lasert & ~isnan(TS.Data(inact_sessions(n)).rts)));
        Identitifier{1} = strcat(num2str(TS.Rat),'-',TS.Date{inact_sessions(n)});
        
        Omissions(1,1) = sum(TS.Data(inact_sessions(n)).omitted(~TS.Data(inact_sessions(n)).lasert));
        Omissions(1,2) = sum(TS.Data(inact_sessions(n)).omitted(TS.Data(inact_sessions(n)).lasert));
        Aborted(1,1) = sum(TS.Data(inact_sessions(n)).aborted(~TS.Data(inact_sessions(n)).lasert));
        Aborted(1,2) = sum(TS.Data(inact_sessions(n)).aborted(TS.Data(inact_sessions(n)).lasert)); 
        Trials(1,1) = sum(sum(TS.Data(inact_sessions(n)).offtype(~TS.Data(inact_sessions(n)).lasert,:)));
        Trials(1,2) = sum(sum(TS.Data(inact_sessions(n)).offtype(TS.Data(inact_sessions(n)).lasert,:)));
        
        
    else
        [inacttemp,Sigma(2,:)] = ModelFitStruct(TS.Data(inact_sessions(n)));
        if flip
            IP(2,:) = -1*inacttemp;
        else
            IP(2,:) = inacttemp;
        end   
        %This gets only reaction times for non-nan trials and for non-laser
        %trials
        RT(2,1) = mean(TS.Data(inact_sessions(n)).rts(~TS.Data(inact_sessions(n)).lasert & ~isnan(TS.Data(inact_sessions(n)).rts)));
        
        %This is for laser
        RT(2,2) = mean(TS.Data(inact_sessions(n)).rts(TS.Data(inact_sessions(n)).lasert & ~isnan(TS.Data(inact_sessions(n)).rts)));
        
        Identitifier{2} = strcat(num2str(TS.Rat),'-',TS.Date{inact_sessions(n)});
        
        Omissions(2,1) = sum(TS.Data(inact_sessions(n)).omitted(~TS.Data(inact_sessions(n)).lasert));
        Omissions(2,2) = sum(TS.Data(inact_sessions(n)).omitted(TS.Data(inact_sessions(n)).lasert));
        Aborted(2,1) = sum(TS.Data(inact_sessions(n)).aborted(~TS.Data(inact_sessions(n)).lasert));
        Aborted(2,2) = sum(TS.Data(inact_sessions(n)).aborted(TS.Data(inact_sessions(n)).lasert));
        Trials(2,1) = sum(sum(TS.Data(inact_sessions(n)).offtype(~TS.Data(inact_sessions(n)).lasert,:)));
        Trials(2,2) = sum(sum(TS.Data(inact_sessions(n)).offtype(TS.Data(inact_sessions(n)).lasert,:)));
    end
end

end

    function [IP, Sigma] = Trans_Data(TS,trans_sessions)
    %To determine the pellet ordering A,B,C, we take the
    %indifference points of all the sessions and laser trials for the
    %experimental days (2 and 3). The matrices are days by row and pairs by
    %column
    %trans_sessions is 3 x 3 matrix of sessions: pairs by days
    
    %This combines the sessions to get a single IP
    CableType = cell(3,3);
    for m = 1:3

        %First determine which offers were used in each session. If the
        %offers were different, combine the sessions
        for n = 2:3
            
            Sess_Offers{n - 1} = TS.Data(trans_sessions(m,n)).Offers;
     
            
        end
        
        if any(Sess_Offers{1}(:) ~= Sess_Offers{2}(:))
            display('session did not have matching offers and was skipped')
            IP = nan(3,2,2);
            Sigma = nan(3,2,2);
            return
        else

        end
        
        CableType{m,1} = 'SHAM';
        All_Choice = zeros(size(Sess_Offers{1}));
        for n = 2:3    
 
            All_Choice = All_Choice + choicemat(TS.Data(trans_sessions(m,n)));
            CableType{m,n} = TS.Inact{trans_sessions(m,n)};
            
        end
        
        [IPAll(m), SigmaAll(m)] = ModelFit(All_Choice, Sess_Offers{1});
        
        %This gets the Pellets used
        Pellets(m,:) = TS.Data(trans_sessions(m,n)).Pellets;
        
        
    end
    
    %This finds the new ordering of preferences. c1 is the pair with the
    %the greatest preference difference
    c1 = abs(IPAll) == max(abs(IPAll));
    Pells3 = unique(Pellets);
    flip = false(1,3);
    
    %c is a logical 1 x 3 with the max IP being the only true
if IPAll(c1) < 0 %In this case the preferences are reversed from the baseline expectations
     
    Pellets_New{1} = Pellets{c1,2};
    Pellets_New{3} = Pellets{c1,1};

    Pellets_New{2} = Pells3{~strcmp(Pells3, Pellets_New{1}) & ~strcmp(Pells3, Pellets_New{3})};
    
    flip(c1) = true;
    
        %this finds the new A vs B
    c2 = ~c1 & any(strcmp(Pellets_New{1}, Pellets'),1);
    
    if sum(c2) > 1
        display('error with code: 1')
    end
    
    %If pellet A is the second pellet than flip the IP for 2
    if strcmp(Pellets_New{1}, Pellets(c2,:)) == [true false]
        
    else
         flip(c2) = true;
    end
    
    %the third pair is not either of the first two pairs
    c3 = ~(c1 | c2);
    
    if strcmp(Pellets_New{2}, Pellets(c3,:)) == [true false]
    else
        flip(c3) = true;
    end
else
    
    Pellets_New{1} = Pellets{c1,1};
    Pellets_New{3} = Pellets{c1,2};
    Pells3 = unique(Pellets);
    Pellets_New{2} = Pells3{~strcmp(Pells3, Pellets_New{1}) & ~strcmp(Pells3, Pellets_New{3})};

    c2 = ~c1 & any(strcmp(Pellets_New{1}, Pellets'),1);
    
    if sum(c2) > 1
        display('error with code: 1')
    end
    
    %Here we're looking at the session of A vs B, if pellet A is 
    %the second pellet (old B) than flip the IP for 2
    if strcmp(Pellets_New{1}, Pellets(c2,:)) == [true false]
    else
        flip(c2) = true;
        
    end
    
    %the third pair is not either of the first two pairs
    c3 = ~(c1 | c2);
    
    if strcmp(Pellets_New{2}, Pellets(c3,:))  == [true false]
    else
        flip(c3) = true;
    end 

end

pair_reorder = [find(c1) find(c2) find(c3)];

%the ip matrix is organized as rows are pairs and columns are sham/opto
%This loops through the pair reorders and the days to input the IP matrix
%The final output for transitivity is a 3D matrix: pair x sham/opto x laser off/on
pair = 0;
for o = pair_reorder
    pair = pair + 1;
    for q = 1:3
        hold on
        if any(strcmpi(CableType{o,q},'sham'))
            [IPtemp, Sigma(pair,1,:)] = ModelFitStruct(TS.Data(trans_sessions(o,q)));
            if flip(o)
                IP(pair,1,:) = -1*IPtemp;
                
            else
                IP(pair,1,:) = IPtemp;
            end
        elseif any(strcmpi(CableType{o,q},'opto'))
            [IPtemp, Sigma(pair,2,:)] = ModelFitStruct(TS.Data(trans_sessions(o,q)));
            if flip(o)
                IP(pair,2,:) = -1*IPtemp;
            else
                IP(pair,2,:) = IPtemp;
            end
        else
            display('Cable type is not sham or opto for trans part of code.')
        end    
    end            
end            
 


end

function [IPhat, Sigmahat] = ModelFitStruct(MAT)
 
     use = 2:size(MAT.offtype,2) - 1;
    
     %no laser trial data:
     ChoiceMAT{1} = ([sum(MAT.offtype(~MAT.lasert & MAT.offer,:)); sum(MAT.offtype(~MAT.lasert & ~MAT.offer,:))]);
     %laser trial data:
     ChoiceMAT{2} = ([sum(MAT.offtype(MAT.lasert & MAT.offer,:)); sum(MAT.offtype(MAT.lasert & ~MAT.offer,:))]);
     
     ratio = MAT.Offers(1,:)./MAT.Offers(2,:);
     
     x = -1*log(ratio(use)); 
     
     for k = 1:2
         
        Choice(1,:) = ChoiceMAT{k}(2,use);
        Choice(2,:) = sum(ChoiceMAT{k}(:, use));
   
        [b, ~,~] = glmfit(x, Choice', 'binomial', 'link', 'probit');

        IPhat(k) = (-1*(b(1)/b(2)));
        Sigmahat(k) = (-1/b(2));

     end
         
end

function [IPhat, Sigmahat] = ModelFit(ChoiceMAT, Offers)
    
     use = 2:size(Offers,2) - 1;
    
     ratio = Offers(1,:)./Offers(2,:);
     
     x = -1*log(ratio(use)); 
     
        Choice(1,:) = ChoiceMAT(2,use);
        Choice(2,:) = sum(ChoiceMAT(:, use));
   
        [b, ~,~] = glmfit(x, Choice', 'binomial', 'link', 'probit');

        IPhat = (-1*(b(1)/b(2)));
        Sigmahat = (1/b(2));

end
         




