cpdef double calc_f1(double recall, double precision):
    """
    Static method to calculate the F1 statistic given precision
    and recall (order is unimportant). Definition:
    F1 = (2 * precision * recall) / (precision + recall)
    """
    if max(precision, recall) == 0:
        return 0
    else:
        product = 2 * precision * recall
        summa = precision + recall
        return product / summa
