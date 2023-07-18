TITLE K-DR channel
: from Klee Ficker and Heinemann
: modified to account for Dax et al.
: M.Migliore 1997

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v (mV)
        ek (mV)		: must be explicitely def. in hoc
	celsius		(degC)
	gkdrbar=.003 (mho/cm2)
        vhalfn=13   (mV)
        a0n=0.02      (/ms)
        zetan=-3    (1)
        gmn=0.7  (1)
	nmax=2  (1)
	q10=1
}


NEURON {
	SUFFIX kdr
	USEION k READ ek WRITE ik
        RANGE gkdr,gkdrbar
	GLOBAL ninf,taun
}

STATE {
	n
}

ASSIGNED {
	ik (mA/cm2)
        ninf
        gkdr
        taun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr = gkdrbar*n
	ik = gkdr*(v-ek)

}

INITIAL {
	rates(v)
	n=ninf
}


FUNCTION alpn(v(mV)) {
  alpn = exp(1.e-3*zetan*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betn(v(mV)) {
  betn = exp(1.e-3*zetan*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))) 
}

DERIVATIVE states {     : exact when v held constant; integrates over dt step
        rates(v)
        n' = (ninf - n)/taun
}

PROCEDURE rates(v (mV)) { :callable from hoc
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1+a)
        taun = betn(v)/(qt*a0n*(1+a))
	if (taun<nmax) {taun=nmax}
}














