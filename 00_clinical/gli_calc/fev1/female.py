
def median(np, age, height, spline):
    return np.exp(
        -10.901689
        + 2.385928 * np.log(height)
        - 0.076386 * np.log(age)
        + spline
    )
