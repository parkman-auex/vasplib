function [HcDesignData,H,T,Hc,time] = PpmsHcData_gen(PPMS_HC_DATA,n)
    T = PPMS_HC_DATA.Temperature_K;
    Hc = PPMS_HC_DATA.SampleHC_muJ_d_K;
    Hc = Hc*1e-6/n;
    H = PPMS_HC_DATA.Magnetic_Field_Oe;
    time = PPMS_HC_DATA.Time_Stamp_sec;
    HcDesignData  = table(time,T,Hc,H);
end