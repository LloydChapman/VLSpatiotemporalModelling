function acc_rate=GetAccRates(rslts)
load(rslts);
acc_rate_A=(acc_Arem+acc_Amov+acc_Aadd)/(acc_Arem+acc_Amov+rej_Aadd+rej_Arem+rej_Amov+acc_Aadd);
acc_rate=[acc_rate_p,acc_rate_E,acc_rate_I,acc_rate_ERmove,acc_rate_R,acc_rate_AIRmove,acc_rate_AR,acc_rate_RLNO,acc_rate_RLO,acc_rate_A,acc_rate_PA];
