function [spl_shift] = checkIfSPLShift(rabbit, session)
% Function that determines if the session is during the time period when
% there was an error in our stimulus presentation that presented the
% stimuli 3dB higher than it was supposed to. This was fixed on 8/5/22
% J. Fritzinger, 9/21/22

spl_shift = 0;
switch rabbit
	case 22
		spl_shift = 3;
	case 23
		spl_shift = 3;
    case 24
        spl_shift = 3;
    case 25
        if session < 676 
            spl_shift = 3;
        end
    case 27
        if session < 81
            spl_shift = 3;
        end
end


end