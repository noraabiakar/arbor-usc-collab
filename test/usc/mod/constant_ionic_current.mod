TITLE constant ionic current

NEURON {
    SUFFIX constant_ionic_current
    USEION x READ ex WRITE ix
}

UNITS {
    (mA) = (milliamp)
}

PARAMETER {
    current = 0 (mA)
}


STATE {
}

ASSIGNED {
}

INITIAL {
}

BREAKPOINT {
    ix = current
}
