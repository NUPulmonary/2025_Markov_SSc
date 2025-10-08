
def median(np, age, height, spline):
    return np.exp(
        -10.5861
        + 0.1433 * np.log(age)
        + 2.3155 * np.log(height)
        + spline
    )
