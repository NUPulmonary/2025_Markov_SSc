
def median(np, age, height, spline):
    return np.exp(
        -11.399108
        + 2.462664 * np.log(height)
        - 0.011394 * np.log(age)
        + spline
    )
