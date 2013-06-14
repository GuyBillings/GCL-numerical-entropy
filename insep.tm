double InSepC P((int,int,int,int,int,double,double,double,int,int));

:Begin:
:Function:       insepc
:Pattern:        InSepC[tc_Integer, mu_Integer, dini_Integer, dinc_Integer, dmax_Integer, pini_Real, pinc_Real, pmax_Real, gamma_Integer, prec_Integer]
:Arguments:      { tc, mu, dini, dinc, dmax, pini, pinc, pmax, gamma, prec }
:ArgumentTypes:  { Integer, Integer, Integer, Integer, Integer, Real, Real, Real, Integer, Integer }
:ReturnType:     Manual
:End:

:Evaluate: InSepC::usage = "InfoPortraitC[tc_,m_Integer,dini_Integer,dinc_Integer,dmax_Integer,pini_Real,pinc_Real,pmax_Real,g_Integer,bitprecision_Integer]"