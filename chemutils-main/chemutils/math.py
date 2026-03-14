import numpy as np
import numpy.typing as npt


def get_signed_dihedral_angle(
    a: npt.NDArray[np.float32],
    b: npt.NDArray[np.float32],
    c: npt.NDArray[np.float32],
    d: npt.NDArray[np.float32],
    /,
) -> npt.NDArray[np.float32]:
    """
    Calculate the signed dihedral angle in radians formed by the points a, b, c and d.

    This is the signed angle between planes formed by the points a, b, c and b, c, d, as viewed down the axis
    from b to c.

    Args:
        a: Coordinates of the first point
        b: Coordinates of the second point
        c: Coordinates of the third point
        d: Coordinates of the fourth point

    Returns:
        The dihedral angle in radians.
    """
    ab = b - a
    bc = c - b  # Vector from point b to c
    cd = d - c  # Vector from point c to d

    # Normal vector of plane formed by a, b, c
    n1 = np.cross(ab, bc)
    # Normal vector of plane formed by b, c, d
    n2 = np.cross(bc, cd)

    atanarg_1 = np.linalg.norm(bc) * np.dot(ab, np.cross(bc, cd))

    atanarg_2 = np.dot(n1, n2)

    # TODO not sure why np.arctan2 does not have type annotations

    angles: npt.NDArray[np.float32] = np.arctan2(atanarg_1, atanarg_2)
    return angles
