"""Meshed rigid punches for CalculiX (mm): a cylinder ring or a hemispherical
shell of C3D8 elements, whose every node is prescribed (so the punch is
kinematically rigid) and whose outer face is the master contact surface."""
import numpy as np

FACES = [(0, 1, 2, 3), (4, 7, 6, 5), (0, 4, 5, 1), (1, 5, 6, 2), (2, 6, 7, 3), (3, 7, 4, 0)]


def hex_volume_sign(x):
    """Sign of the Jacobian at the centre of a trilinear hex (8 x 3)."""
    d1 = (x[1] + x[2] + x[5] + x[6] - x[0] - x[3] - x[4] - x[7])
    d2 = (x[2] + x[3] + x[6] + x[7] - x[0] - x[1] - x[4] - x[5])
    d3 = (x[4] + x[5] + x[6] + x[7] - x[0] - x[1] - x[2] - x[3])
    return np.sign(np.dot(d1, np.cross(d2, d3)))


def _finish(nodes, conn, outer_test):
    nodes = np.asarray(nodes, float)
    conn = np.asarray(conn, int)
    for e in range(len(conn)):
        if hex_volume_sign(nodes[conn[e]]) < 0:
            conn[e] = conn[e][[4, 5, 6, 7, 0, 1, 2, 3]]
    faces = []
    for e in range(len(conn)):
        for f, loc in enumerate(FACES):
            if all(outer_test(nodes[conn[e][i]]) for i in loc):
                faces.append((e, f + 1))
                break
    return nodes, conn, faces


def cylinder_ring(center, R, y0, y1, theta_max_deg=100.0, n_theta=200, thick=0.5):
    """Ring of the cylinder with axis along y through center (x, z): angles from
    the downward vertical in [-theta_max, theta_max]; outer face r = R."""
    th = np.radians(np.linspace(-theta_max_deg, theta_max_deg, n_theta + 1))
    ys = [y0, y1]
    rs = [R, R - thick]
    idx = {}
    nodes = []
    for k, r in enumerate(rs):
        for j, y in enumerate(ys):
            for i, t in enumerate(th):
                idx[i, j, k] = len(nodes)
                nodes.append([center[0] + r * np.sin(t), y, center[1] - r * np.cos(t)])
    conn = []
    for i in range(n_theta):
        conn.append([idx[i, 0, 0], idx[i + 1, 0, 0], idx[i + 1, 1, 0], idx[i, 1, 0],
                     idx[i, 0, 1], idx[i + 1, 0, 1], idx[i + 1, 1, 1], idx[i, 1, 1]])
    c = np.array([center[0], center[1]])
    return _finish(nodes, conn,
                   lambda p: abs(np.hypot(p[0] - c[0], p[2] - c[1]) - R) < 1e-6 * R)


def sphere_shell(center, R, n=24, thick=0.5):
    """Lower hemisphere of the sphere of radius R about center (cubed-sphere
    mapping: the bottom face n x n and the lower halves of the four sides
    n x n/2), extruded inward by thick; outer face r = R."""
    pts = {}
    nodes = []

    def node(p):
        key = tuple(np.round(p, 9))
        if key not in pts:
            pts[key] = len(nodes)
            nodes.append(p)
        return pts[key]

    u = np.linspace(-1, 1, n + 1)
    half = np.linspace(-1, 0, n // 2 + 1)
    patches = []
    # bottom face z = -1
    patches.append([[np.array([a, b, -1.0]) for a in u] for b in u])
    # sides x = +-1, y = +-1, z from -1 to 0
    for s in (-1.0, 1.0):
        patches.append([[np.array([s, a, c]) for a in u] for c in half])
        patches.append([[np.array([a, s, c]) for a in u] for c in half])
    conn = []
    c0 = np.asarray(center, float)
    for P in patches:
        m, l = len(P), len(P[0])
        ids = {}
        for a in range(m):
            for b in range(l):
                d = P[a][b] / np.linalg.norm(P[a][b])
                ids[a, b, 0] = node(c0 + R * d)
                ids[a, b, 1] = node(c0 + (R - thick) * d)
        for a in range(m - 1):
            for b in range(l - 1):
                conn.append([ids[a, b, 0], ids[a, b + 1, 0], ids[a + 1, b + 1, 0], ids[a + 1, b, 0],
                             ids[a, b, 1], ids[a, b + 1, 1], ids[a + 1, b + 1, 1], ids[a + 1, b, 1]])
    return _finish(nodes, conn, lambda p: abs(np.linalg.norm(p - c0) - R) < 1e-6 * R)
