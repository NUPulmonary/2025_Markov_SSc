
def median(np, age, height, spline):
    return np.exp(
        0.9189568
        - 0.1840671 * np.log(height)
        - 0.0461306 * np.log(age)
        + spline
    )
