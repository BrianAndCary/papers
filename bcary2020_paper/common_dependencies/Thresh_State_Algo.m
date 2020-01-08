function [output_states] = Thresh_State_Algo(state_times,features,labels,thresholds)

% features should be matrix with each col corrresponding to vect of feature
% data. The order of features should be:
% dratio, tratio, emg, movt

% thresholds are in same order as feature data i.e.:

t = state_times;
d_to_b_ratio = features(:,1);
t_to_d_ratio = features(:,2);
emg_average = features(:,3);
move_average = features(:,4);

nrem_thresh = thresholds(:,1);
rem_thresh = thresholds(:,2);
emg_thresh = thresholds(:,3);
movt_thresh = thresholds(:,4);

% brain state classifier
% 1 = awake
% 2 = rem
% 3 = nrem

brain_state = labels;
if isempty(brain_state)
    
    % first classification point
    if (move_average(1) > movt_thresh*1.2) || (emg_average(1) > emg_thresh*1.2)
         brain_state(1) = 1;
    elseif t_to_d_ratio(1) > rem_thresh*2
        brain_state(1) = 2;
    elseif d_to_b_ratio(1) > nrem_thresh
        brain_state(1) = 3;
    else
        brain_state(1) = 1;
    end
    
    % second classification point
    if (move_average(2) > movt_thresh*1.2) || (emg_average(1) > emg_thresh*1.2)
         brain_state(2) = 1;
    elseif t_to_d_ratio(2) > rem_thresh*2
        brain_state(2) = 2;
    elseif d_to_b_ratio(2) > nrem_thresh
        brain_state(2) = 3;
    else
        brain_state(2) = 1;
    end
    
    start_loop = 3;
else
    start_loop = 1;
end

% classification loop
for min = start_loop:length(t)
    one_back = brain_state(end);
    two_back = brain_state(end-1);
    if two_back == 1
        if one_back == 1
            if (move_average(min) > movt_thresh*0.8) || (emg_average(1) > emg_thresh*0.8)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*3
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*1.1
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > movt_thresh*1.1) || (emg_average(1) > emg_thresh*1.1)
                brain_state(end) = 1;
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*1.1
                try
                if (move_average(min-1) > movt_thresh*0.25) || (emg_average(1) > emg_thresh*0.25)
                    brain_state(end) = 1;
                end                    
                brain_state(end+1) = 3;
                end
            else
                brain_state(end) = 1;
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > movt_thresh*1.2) || (emg_average(1) > emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*1.1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end
    
    if two_back == 2
        if one_back == 1
            if (move_average(min) > movt_thresh*1) || (emg_average(1) > emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*1
                try
                if (move_average(min-1) < movt_thresh*0.9) || (emg_average(1) < emg_thresh*1.2)
                    brain_state(end) = 2;
                end
                brain_state(end+1) = 2;
                end
            elseif d_to_b_ratio(min) > nrem_thresh*1            
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > movt_thresh*1.2) || (emg_average(1) > emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*0.85
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*0.85
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > movt_thresh*1.1) || (emg_average(1) > emg_thresh*1.1)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*1.1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end     
        

    if two_back == 3
        if one_back == 1
            if (move_average(min) > movt_thresh*1) || (emg_average(1) > emg_thresh*1)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*1
                rem_ratio = t_to_d_ratio(min-1)/rem_thresh;
                nrem_ratio = d_to_b_ratio(min-1)/nrem_thresh;
                
                try
                if (move_average(min-1) < movt_thresh*0.9) || (emg_average(1) < emg_thresh*0.9)
                    if rem_ratio > nrem_ratio
                        brain_state(end) = 2;
                    else
                        brain_state(end) = 3;
                    end
                end
                catch
                    if rem_ratio > nrem_ratio
                        brain_state(end) = 2;
                    else
                        brain_state(end) = 3;
                    end
                end
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*1
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end 
        
        if one_back == 2
            if (move_average(min) > movt_thresh*1.2) || (emg_average(1) > emg_thresh*1.2)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*0.9
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*0.85
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end         
        
        if one_back == 3
            if (move_average(min) > movt_thresh*1.25) || (emg_average(1) > emg_thresh*1.25)
                brain_state(end+1) = 1;
            elseif t_to_d_ratio(min) > rem_thresh*1
                brain_state(end+1) = 2;
            elseif d_to_b_ratio(min) > nrem_thresh*0.9
                brain_state(end+1) = 3;
            else
                brain_state(end+1) = 1;
            end
        end                 
    end
    
end

output_states = brain_state(end-length(t)+1:end);

end