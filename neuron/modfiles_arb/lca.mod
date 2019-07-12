TITLE l-calcium channel
: l-type calcium channel


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
    :KTOMV = .0853 (mV/degC)
}

PARAMETER {
    v (mV)
    celsius 	(degC)
    glcabar		 (mho/cm2)
    ki=.001 (mM)
    tfa=1

    cai : NEURON
    cao : NEURON
}


NEURON {
    THREADSAFE
    SUFFIX lca
    USEION lca READ elca WRITE ilca VALENCE 2
    USEION ca READ cai, cao VALENCE 2
    RANGE glcabar, minf, matu

    RANGE elca, ilca, cai : NEURON
}

STATE {
    m
}

ASSIGNED {
    glca (mho/cm2)
    minf
    matu   (ms)

    f   : NEURON
    nu  : NEURON
    ghk : NEURON

    elca : NEURON
    ilca : NEURON
}

INITIAL {
    rate(v)
    m = minf
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    glca = glcabar*m*m*h2(cai)

    f = KTF(celsius)/2
    nu = v/f
    ghk=-f*(1. - (cai/cao)*exp(nu))*efun(nu) : NEURON
    ilca = glca*ghk

}

FUNCTION h2(cai) {
    h2 = ki/(ki+cai)
}

FUNCTION KTF(celsius) {
    KTF = ((25./293.15)*(celsius + 273.15))
}

FUNCTION alp(v) {
    alp = 15.69*(-1.0*v+81.5)/(exp((-1.0*v+81.5)/10.0)-1.0)
}

FUNCTION bet(v) {
    bet = 0.29*exp(-v/10.86)
}

DERIVATIVE state {  
    rate(v)
    m' = (minf - m)/matu
}

PROCEDURE rate(v) {
    LOCAL a
    a = alp(v)
    matu = 1/(tfa*(a + bet(v)))
    minf = tfa*a*matu
}
 
FUNCTION efun(z) { : NEURON
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

