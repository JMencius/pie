import math

def safediv(a: float, b: float):
    if (b == 0) or math.isnan(b):
        return math.nan
    else:
        return a / b
