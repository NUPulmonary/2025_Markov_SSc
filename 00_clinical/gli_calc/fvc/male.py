
def median(np, age, height, spline):
    return np.exp(
        -12.629131
        + 2.727421 * np.log(height)
        + 0.009174 * np.log(age)
        + spline
    )
