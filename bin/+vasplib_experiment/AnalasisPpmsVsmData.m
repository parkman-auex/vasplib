function AnalasisPpmsVsmData(PpmsVsmData)
    % time
    time= PpmsVsmData.Time_Stamp_sec;
    time= time-time(1);
    Duration =  time(end);
    time_step = time(6)-time(6);
    [th,tm,ts]=sec2hms(Duration);
    fprint('Your sample has been measured about d% hours d% minutes and %f seconds.\n',th,tm,ts);
    % find data cut

end
function [th,tm,ts]=sec2hms(sec)
    th = floor(sec,3600);
    tm = floor(sec-th*3600,60);
    ts = sec -th*3600-tm*60;
end



