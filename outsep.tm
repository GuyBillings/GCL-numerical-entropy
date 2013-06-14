double OutSepC P((int,int,int,int,int,double,double,double,int,int));

:Begin:
:Function:       outsepc
:Pattern:        OutSepC[tc_Integer, mu_Integer, dini_Integer, dinc_Integer, dmax_Integer, pini_Real, pinc_Real, pmax_Real, gamma_Integer, prec_Integer]
:Arguments:      { tc, mu, dini, dinc, dmax, pini, pinc, pmax, gamma, prec }
:ArgumentTypes:  { Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Integer, Integer }
:ReturnType:     Manual
:End:

:Evaluate: OutSepC::usage = "InfoPortraitC[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer,bitprecision_Integer]"