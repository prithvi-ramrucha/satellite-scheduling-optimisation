def invervals_overlap(a, b):
    """It two time intervals a and overlap, the function returns True.
    Each interval is represented as a tuple (start, end)."""
    return a[1] > b[0] and a[0] < b[1]

def compute_overlaps(windows):
    """Computes the overlapping windows from a list of windows. Each window is represented as a tuple 
    (satellite_name, start_time, end_time, duration, capacity). A list of tuples (i, j) is returned, 
    where i and j are the indices of the overlapping windows."""
    O = []
    for i in range(len(windows)):
        for j in range(i+1, len(windows)):
            a = windows[i][1:3]
            b = windows[j][1:3]
            if invervals_overlap(a, b):
                O.append((i, j))
    return O