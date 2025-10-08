
def median(np, age, height, spline):
    return np.exp(
        -7.034920
        - 0.012425 * np.log(age)
        + 2.018368 * np.log(height)
        + spline
    )
