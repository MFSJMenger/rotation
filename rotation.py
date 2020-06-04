import numpy as np


def generate_translation(p):
    """Generate translation matrix as
        ( 1  0  0  -x )
    T = ( 0  1  0  -y )
        ( 0  0  1  -z )
        ( 0  0  0   1 )

        ( 1  0  0   x )
   TT = ( 0  1  0   y )
        ( 0  0  1   z )
        ( 0  0  0   1 )

    """
    t = np.identity(4)
    tt = np.identity(4)
    for i in range(3):
        t[i][3] = -p[i]
        tt[i][3] = p[i]
    return t, tt


def rotation_matrix_xaxis(U):
    """
    if already in x axis: Rx = Identity

        ( 1   0    0  0 )
   Rx = ( 0  c/d -b/d 0 )
        ( 0  b/d  c/d 0 )
        ( 0   0    0  1 )

        ( 1   0    0  0 )
  Rxt = ( 0  c/d  b/d 0 )
        ( 0 -b/d  c/d 0 )
        ( 0   0    0  1 )

    """
    b = U[1]
    c = U[2]
    d = np.sqrt(b**2 + c**2)
    Rx = np.identity(4)
    Rxt = np.copy(Rx)
    if np.round(d, 6) == 0.0:
        # we are already in the x axis
        return out
    Rx[1:3, 1:3] = np.array([[c/d, -b/d], [b/d, c/d]])
    Rxt[1:3, 1:3] = np.array([[c/d, b/d], [-b/d, c/d]])
    return Rx, Rxt


def rotation_matrix_yaxis(U):
    """
        ( d  0 -a  0 )
   Ry = ( 0  1  0  0 )
        ( a  0  d  0 )
        ( 0  0  0  1 )

        ( d  0  a  0 )
  Ryt = ( 0  1  0  0 )
        (-a  0  d  0 )
        ( 0  0  0  1 )

    """
    a = U[0]
    b = U[1]
    c = U[2]
    d = np.sqrt(b**2 + c**2)
    Ry = np.identity(4)
    Ry[0][0] = d
    Ry[2][2] = d
    Ryt = np.copy(Ry)
    Ry[0][2] = -a
    Ry[2][0] = a
    Ryt[0][2] = a
    Ryt[2][0] = -a
    return Ry, Ryt

def rotation_matrix_zaxis(theta):
    """
        ( cos(t) -sin(t)  0  0 )
   Ry = ( sin(t)  cos(t)  0  0 )
        ( 0        0      0  0 )
        ( 0        0      0  1 )

    """
    out = np.identity(4)
    out[:2,:2] = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta),  np.cos(theta)]])
    return out

def generate_rotation_operator(p1, p2, theta):
    """
        Generate rotation around arbitray axis spanned by two points
        p1, p2:

        R = Tt * Rxt * Ryt * Rz * Ry * Rx * T
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    p = norm_difference(p1, p2)

    T, Ttrans = generate_translation(p1)
    Rx, Rxt = rotation_matrix_xaxis(p)
    Ry, Ryt = rotation_matrix_yaxis(p)
    Rz = rotation_matrix_zaxis(theta)
    return Ttrans.dot(Rxt.dot(Ryt.dot(Rz.dot(Ry.dot(Rx.dot(T))))))

def translate(p, T):
    return T.dot(p)

def norm(vec):
    return np.sqrt(vec.dot(vec))
