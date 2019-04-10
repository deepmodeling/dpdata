import numpy as np


def posi_diff(box, r0, r1) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    p0 = (np.dot(rbox, r0))
    p1 = (np.dot(rbox, r1))
    dp = p0 - p1
    shift = np.zeros(3)
    for dd in range(3) :
        if dp[dd] >= 0.5 : 
            dp[dd] -= 1
        elif dp[dd] < -0.5 :
            dp[dd] += 1
    dr = np.dot(box.T, dp)    
    return dr


def posi_shift(box, r0, r1) :
    rbox = np.linalg.inv(box)
    rbox = rbox.T
    p0 = (np.dot(rbox, r0))
    p1 = (np.dot(rbox, r1))
    dp = p0 - p1
    shift = np.zeros(3)
    for dd in range(3) :
        if dp[dd] >= 0.5 : 
            shift[dd] -= 1
        elif dp[dd] < -0.5 :
            shift[dd] += 1
    return shift

