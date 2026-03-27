LOGN = 9 # HAWK parameter 8 = HAWK-256, 9 = HAWK-512, 10 = HAWK-1024


MODE = "ORIGIN" # Switch between ORIGIN, PURE, LITE, ONE

# OWL:ORIGIN
ORIGIN = {
    "LAMBDA": 128,
    "C": 7,
    "K": 22,
    "ROUND": 84
}


PURE = {
    "LAMBDA": 128,
    "C": 15,
    "ROUND": 64
}

# OWL:LITE
LITE = {
    "LAMBDA": 128,
}

# OWL:ONE
ONE = {
    "LAMBDA": 128,
}

def OWL_PARAMS(MODE):
    if MODE == "ORIGIN":
        return ORIGIN
    elif MODE == "PURE":
        return PURE
    elif MODE == "LITE":
        return LITE
    elif MODE == "ONE":
        return ONE
    else:
        raise Exception("Unrecognized MODE.")