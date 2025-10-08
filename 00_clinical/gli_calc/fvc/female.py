
def median(np, age, height, spline):
    return np.exp(
        -12.055901
        + 2.621579 * np.log(height)
        - 0.035975 * np.log(age)
        + spline
    )
