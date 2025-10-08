
def median(np, age, height, spline):
    return np.exp(
        -5.159451
        - 0.015390 * np.log(age)
        + 1.618697 * np.log(height)
        + spline
    )
