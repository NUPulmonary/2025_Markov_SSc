
def median(np, age, height, spline):
    return np.exp(
        -10.1128
        + 0.1062 * np.log(age)
        + 2.2259 * np.log(height)
        + spline
    )
