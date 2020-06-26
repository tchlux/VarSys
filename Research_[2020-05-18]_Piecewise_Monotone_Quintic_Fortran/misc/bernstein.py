def bernstein_from_derivatives(xa, xb, ya, yb):
    ya, yb = np.asarray(ya), np.asarray(yb)
    if ya.shape[1:] != yb.shape[1:]:
        raise ValueError('ya and yb have incompatible dimensions.')

    dta, dtb = ya.dtype, yb.dtype
    dt = np.float_

    na, nb = len(ya), len(yb)
    n = na + nb

    c = np.empty((na+nb,) + ya.shape[1:], dtype=dt)

    # compute coefficients of a polynomial degree na+nb-1
    # walk left-to-right
    for q in range(0, na):
        c[q] = ya[q] / spec.poch(n - q, q) * (xb - xa)**q
        for j in range(0, q):
            c[q] -= (-1)**(j+q) * comb(q, j) * c[j]

    # now walk right-to-left
    for q in range(0, nb):
        c[-q-1] = yb[q] / spec.poch(n - q, q) * (-1)**q * (xb - xa)**q
        for j in range(0, q):
            c[-q-1] -= (-1)**(j+1) * comb(q, j+1) * c[-q+j]

    return c
