
def median(np, age, height, spline):
    return np.exp(
        1.022608
        - 0.218592 * np.log(height)
        - 0.027584 * np.log(age)
        + spline
    )
