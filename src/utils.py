
import numpy as np
from itertools import compress

def meanAbsDev(dat):
    return np.mean(np.absolute(dat - np.mean(dat)))

def medianAbsDev(dat):
    median_x = np.median(dat)
    return np.median([np.abs(dat - median_x) for x in dat])


# return list of bools of vals in list which are not outliers
def outliers_modified_z_score(ys):
    threshold = 3.5

    median_y = np.median(ys)
    median_absolute_deviation_y = medianAbsDev(ys)

    if median_absolute_deviation_y == 0:
        return [True] * len(ys)

    modified_z_scores = [0.6745 * (y - median_y) / median_absolute_deviation_y
                         for y in ys]

    ret = [val <= threshold for val in modified_z_scores]

    return ret

# return list of bools of vals in list which are not outliers
def outliers_iqr(ys):
    quartile_1 = np.percentile(ys, 25)
    quartile_3 = np.percentile(ys, 75)

    iqr = quartile_3 - quartile_1
    lower_bound = quartile_1 - (iqr * 1.5)
    upper_bound = quartile_3 + (iqr * 1.5)

    ret = [(lower_bound <= val <= upper_bound) for val in ys]

    return ret

# return all non outliers in list based off user specified outlier test method
def rmOutliers(dat, _outlierTest = 1):
    if _outlierTest == 0:
        return dat
    elif _outlierTest == 1:
        booList = outliers_iqr(dat)
        return list(compress(dat, booList)), booList
    elif _outlierTest == 2:
        booList = outliers_modified_z_score(dat)
        return list(compress(dat, booList)), booList
    else:
        raise Exception(str(_outlierTest) + " is not a valid option")

