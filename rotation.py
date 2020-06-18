"""Rotate a point about an arbitrary axis (3 dimensions)

according to: http://paulbourke.net/geometry/rotate/
"""
from copy import deepcopy
#
import numpy as np


def generate_rotation_operator(p1, p2, theta):
    """Generate rotation function around an arbitray axis spanned by two points p1, p2:

        1) translate space so that the rotation axis passes through the origion
           by translating the system so that `p1` is the origin using the
           translation Matrix `T`
        
        2) Rotate space about the x-axis so that the rotation axis lies in the xz plane
           using the transformation matrix `Rx`

        3) rotate space about the y axis so that the rotation axis lies along the z-axis
           using the transformation matrix `Ry`

        4) perform the desired rotation by theta about the z-axis
           using rotation matrix `Rz`

        5) apply the inverse of step (3) using the inverse of `Ry`

        6) apply the inverse of step (2) using the inverse of `Rx`

        7) apply the inverse of step (1) using the inverse of `T`

    The full rotation matrix is then given as:

        R = Tt * Rxt * Ryt * Rz * Ry * Rx * T

    Parameters
    ----------
    p1, p2: np.array, dim(3)
        points that define the axis

    theta: float, degree
        rotation angle in degree

    Returns
    -------
    function
        rotation function for this particular rotation
    """
    # from degree
    theta = theta/180.0*np.pi
    # 
    p1, p2 = np.array(p1), np.array(p2)
    p = norm_difference(p1, p2)
    #
    T, Ttrans = generate_translation(p1)
    Rx, Rxt = rotation_matrix_xaxis(p)
    Ry, Ryt = rotation_matrix_yaxis(p)
    Rz = rotation_matrix_zaxis(theta)
    # 
    return transformation_operation(Ttrans.dot(Rxt.dot(Ryt.dot(Rz.dot(Ry.dot(Rx.dot(T)))))))




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

    Parameters
    ----------
    p: np.array(3)
        displacement vector
    
    Returns
    -------
    np.array((4,4))
        Transformation matrix

    """
    t = np.identity(4)
    tt = np.identity(4)
    for i in range(3):
        t[i][3] = -p[i]
        tt[i][3] = p[i]
    return t, tt


def rotation_matrix_xaxis(U):
    """Generate roation matrix to x-axis
    if already in x axis: Rx = Identity

          ( 1   0    0  0 )
     Rx = ( 0  c/d -b/d 0 )
          ( 0  b/d  c/d 0 )
          ( 0   0    0  1 )
   
          ( 1   0    0  0 )
    Rxt = ( 0  c/d  b/d 0 )
          ( 0 -b/d  c/d 0 )
          ( 0   0    0  1 )

    Parameters
    ----------
    U: np.array(3)
    
    Returns
    -------
    np.array((4,4))
        Transformation matrix
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
    """Generate rotation matrix around y axis

          ( d  0 -a  0 )
     Ry = ( 0  1  0  0 )
          ( a  0  d  0 )
          ( 0  0  0  1 )
    
          ( d  0  a  0 )
    Ryt = ( 0  1  0  0 )
          (-a  0  d  0 )
          ( 0  0  0  1 )

    Parameters
    ----------
    U: np.array(3)
    
    Returns
    -------
    np.array((4,4))
        Transformation matrix
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
    """Creates rotation matrix around x axis according to:

        ( cos(t) -sin(t)  0  0 )
   Ry = ( sin(t)  cos(t)  0  0 )
        ( 0        0      0  0 )
        ( 0        0      0  1 )

    Parameters
    ----------
    theta: float, radians
        angle in radians
    
    Returns
    -------
    np.array((4,4))
        Rotation matrix

    """
    out = np.identity(4)
    out[:2,:2] = np.array([[np.cos(theta), -np.sin(theta)],
                           [np.sin(theta),  np.cos(theta)]])
    return out


def norm_difference(p1, p2):
    p = p2 - p1
    return p/np.sqrt(p.dot(p))


def transformation_operation(T):
    """Generate a transformation function from a quaterion matrix
    Parameters
    ----------
    T: np.array((4,4))
        Transformation Matrix

    Returns
    -------
    function
        function to perform transformation in 3D
    """

    code = f"""
def _helper(p):
    out = np.zeros(3)
    #
    out[0] = p[0]*{T[0,0]} + p[1]*{T[0,1]} + p[2]*{T[0,2]} + {T[0, 3]}
    out[1] = p[0]*{T[1,0]} + p[1]*{T[1,1]} + p[2]*{T[1,2]} + {T[1, 3]}
    out[2] = p[0]*{T[2,0]} + p[1]*{T[2,1]} + p[2]*{T[2,2]} + {T[2, 3]}
    #
    return out
    """
    dct = {'np': np}
    exec(code, dct)
    #
    return dct['_helper']
