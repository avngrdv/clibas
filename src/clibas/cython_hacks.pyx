# idx_finder.pyx
import numpy as np

cimport numpy as np


def idx_finder(np.ndarray[np.uint8_t, ndim=1] arr,
               np.ndarray[np.uint8_t, ndim=1] left,
               np.ndarray[np.uint8_t, ndim=1] right,
               int tol):

    cdef int left_size = left.shape[0]
    cdef int right_size = right.shape[0]
    cdef int n = arr.shape[0]
    cdef int i, j, mis
    cdef int idx1 = -1
    cdef int idx2 = -1

    # search forward for idx1 (first match with left)
    for i in range(n - left_size + 1):
        mis = 0
        for j in range(left_size):
            if arr[i + j] != left[j]:
                mis += 1
                if mis > tol:
                    break
        if mis <= tol:
            idx1 = i
            break  # stop as soon as the first match is found

    # search backward for idx2 (last match with right)
    for i in range(n - right_size, -1, -1):
        mis = 0
        for j in range(right_size):
            if arr[i + j] != right[j]:
                mis += 1
                if mis > tol:
                    break
        if mis <= tol:
            idx2 = i
            break  # stop as soon as the last match found scanning backward

    return idx1 if idx1 != -1 else None, (idx2 + right_size) if idx2 != -1 else None

