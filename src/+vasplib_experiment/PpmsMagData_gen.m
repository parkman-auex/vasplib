function [MagDesignData,M,H,T,time] = PpmsMagData_gen(PPMS_VSM_DATA)
    T = PPMS_VSM_DATA.Temperature_K;
    M = PPMS_VSM_DATA.Moment_emu;
    H = PPMS_VSM_DATA.Magnetic_Field_Oe;
    time = PPMS_VSM_DATA.Time_Stamp_sec;
    MagDesignData  = table(time,M,H,T);
end