%This code will relabel the triggers for the Filtered sounds from 10 to 20


for x=1:length(EEG.event)
    EEG.event(x).type=20;
end
