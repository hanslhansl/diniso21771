"""Implementation of DIN ISO 21771 - Gears – Cylindrical involute gears and gear pairs – Concepts and geometry.

See https://github.com/hanslhansl/diniso21771."""

from typing import Optional, Any
from scipy import optimize
import math as m

def involute(alpha : float):
    """alpha in degrees"""
    return m.tan(m.radians(alpha)) - m.radians(alpha)
def inverse_involute(alpha : float, anfangswert = 20.):
    try:
        return float(optimize.newton(lambda x: involute(x) - alpha, anfangswert))
    except RuntimeError:
        assert(False)

Ritzel = 0
Rad = 1
_indices = (Ritzel, Rad)

def d(z : int, m_t : float):
    return abs(z) * m_t
def d_b(d : float, alpha_t : float):
    return d * m.cos(m.radians(alpha_t))
def d_a(m_n : float, z : int, x : float, h_aP : float, k : float, d : float):
    return d + 2 * z / abs(z) * (x * m_n + h_aP + k * m_n)
def d_f(m_n : float, z : int, x : float, h_fP : float, d : float):
    return d - 2 * z / abs(z) * (h_fP - x * m_n)
def d_w(d_b : float, alpha_wt : float):
    return d_b / m.cos(m.radians(alpha_wt))
def h(m_n : float, k : float, h_aP : float, h_fP : float):
    return h_aP + k * m_n + h_fP
def z_min(m_n : float, x : float, h_aP0, alpha_t : float, beta : float):
    return 2 * m.cos(m.radians(beta)) * (h_aP0 / m_n - x) / m.sin(m.radians(alpha_t))**2
def gamma(z : int, x : float, alpha_n : float, alpha_t : float):
    return inverse_involute( m.pi / 2 / z + 2 * x / z * m.tan(m.radians(alpha_n)) + involute(alpha_t))
def d_amax(z : int, m_t : float, alpha_t : float, gamma : float):
    return m_t * z * m.cos(m.radians(alpha_t)) / m.cos(m.radians(gamma))
def alpha_t(alpha_n : float, beta : float):
    return m.degrees(m.atan(m.tan(m.radians(alpha_n)) / m.cos(m.radians(beta))))
def alpha_wt(z : tuple[int, int], x : tuple[float, float], alpha_n : float, alpha_t : float):
    return inverse_involute(involute(alpha_t) + 2 * sum(x) / sum(z) * m.tan(m.radians(alpha_n)) )
def u(z : tuple[int, int]):
    return z[Rad] / z[Ritzel]
def m_t(m_n : float, beta : float):
    return m_n / m.cos(m.radians(beta))
def beta_b(beta : float, alpha_t : float):
    return m.degrees(m.atan(m.tan(m.radians(beta)) * m.cos(m.radians(alpha_t))))
def a_w(z_2 : int, d_w : float):
    return (d_w[Rad] + z_2 / abs(z_2) * d_w[Ritzel]) / 2
def p_n(m_n : float):
    return m_n * m.pi
def p_t(p_n : float, beta : float):
    return p_n / m.cos(m.radians(beta))
def epsilon_alpha(z : tuple[int, int], d_a : tuple[float, float], d_b : tuple[float, float], p_t : float, alpha_t : float, alpha_wt : float):
    return (m.sqrt(d_a[Ritzel]**2 - d_b[Ritzel]**2) + z[Rad] / abs(z[Rad]) * m.sqrt(d_a[Rad]**2 - d_b[Rad]**2)
                              - (d_b[Ritzel] + z[Rad] / abs(z[Rad]) * d_b[Rad]) * m.tan(m.radians(alpha_wt))) / (2 * p_t * m.cos(m.radians(alpha_t)))
def epsilon_beta(b : float, p_n : float, beta : float):
    return b * m.sin(m.radians(beta)) / p_n
def epsilon_gamma(epsilon_alpha : float, epsilon_beta : float):
    return epsilon_alpha + epsilon_beta

class GearGeometry:
    def __init__(self,
            m_n : float,
            z: tuple[int, int],
            x: tuple[float, float],
            bezugsprofil : tuple[Any, Any],
            beta : float,
            k : int,
            b : Optional[float] = None,
            b_d_1_verhältnis : Optional[float] = None,
            _print = print):
        """
        Parameters:
        - m_n: Modul
        - z: Zähnezahlen. z1 ist immer positiv. Für Außenradpaare ist z2 positiv, für Innenradpaare ist z2 negativ.
        - x: Profilverschiebungsfaktoren
        """

        assert all(z > 0 for z in z), "Innenverzahnung nicht implementiert"

        self.m_n = m_n
        self.z = z
        self.x = x
        self.beta = beta
        self.k = k

        _print("Getriebegeometrie")
        [_print(key, "=", value) for key, value in vars(self).items()]

        assert all(bezugsprofil[idx].alpha_n == bezugsprofil[Ritzel].alpha_n for idx in _indices)
        self.alpha_n = bezugsprofil[Ritzel].alpha_n
        self.h_aP = tuple(bezugsprofil[idx].h_aP_s * self.m_n for idx in _indices)
        self.h_fP = tuple(bezugsprofil[idx].h_fP_s * self.m_n for idx in _indices)
        self.rho_fP = tuple(bezugsprofil[idx].rho_fP_s * self.m_n for idx in _indices)
        self.s_pr = tuple(bezugsprofil[idx].s_pr_s * self.m_n for idx in _indices)
        _print("α_n =", self.alpha_n)
        _print("h_aP =", self.h_aP)
        _print("h_fP =", self.h_fP)
        _print("ρ_fP =", self.rho_fP)

        self.alpha_t = alpha_t(self.alpha_n, self.beta)
        _print("α_t =", self.alpha_t)

        self.alpha_wt = alpha_wt(self.z, self.x, self.alpha_n, self.alpha_t)
        _print("α_wt =", self.alpha_wt)

        self.u = u(self.z)
        _print("u =", self.u)

        self.m_t = m_t(self.m_n, self.beta)
        _print("m_t =", self.m_t)

        self.d = tuple(d(self.z[idx], self.m_t) for idx in _indices)
        _print("d =", self.d)

        self.d_b = tuple(d_b(self.d[idx], self.alpha_t) for idx in _indices)
        _print("d_b =", self.d_b)

        self.d_a = tuple(d_a(self.m_n, self.z[idx], self.x[idx], self.h_aP[idx], self.k, self.d[idx]) for idx in _indices)
        _print("d_a =", self.d_a)

        self.d_f = tuple(d_f(self.m_n, self.z[idx], self.x[idx], self.h_fP[idx], self.d[idx]) for idx in _indices)
        _print("d_f =", self.d_f)

        self.d_w = tuple(d_w(self.d_b[idx], self.alpha_wt) for idx in _indices)
        _print("d_w =", self.d_w)

        if b != None:
            self.b = b
        elif b_d_1_verhältnis != None:
            self.b = self.d[Ritzel] * b_d_1_verhältnis
        else:
            raise ValueError("either b or b_d_1_verhältnis must be specified as argument")
        _print("b =", self.b)

        self.beta_b = beta_b(self.beta, self.alpha_t)
        _print("β_b =", self.beta_b)
        
        self.h = tuple(h(self.m_n, self.k, self.h_aP[idx], self.h_fP[idx]) for idx in _indices)
        _print("h =", self.h)

        self.a_w = a_w(self.z[Rad], self.d_w)
        _print("a_w =", self.a_w)

        # Profilüberdeckung

        self.p_n = p_n(self.m_n)
        _print("p_n =", self.p_n)

        self.p_t = p_t(self.p_n, self.beta)
        _print("p_t =", self.p_t)

        self.epsilon_alpha = epsilon_alpha(self.z, self.d_a, self.d_b, self.p_t, self.alpha_t, self.alpha_wt)
        _print("ε_α =", self.epsilon_alpha)

        self.epsilon_beta = epsilon_beta(self.b, self.p_n, self.beta)
        _print("ε_β =", self.epsilon_beta)

        self.epsilon_gamma = epsilon_gamma(self.epsilon_alpha, self.epsilon_beta)
        _print("ε_γ =", self.epsilon_gamma)
        assert self.epsilon_gamma > 1

        # Unterschnitt

        self.h_aP0 = self.h_fP
        _print("h_aP0 =", self.h_aP0)

        _z_min = tuple(z_min(self.m_n, self.x[idx], self.h_aP0[idx], self.alpha_t, self.beta) for idx in _indices)
        _print("z_min =", _z_min)
        assert self.z[Ritzel] > _z_min[Ritzel]
        assert self.z[Rad] > _z_min[Rad]

        # Spitzwerden

        _gamma = tuple(gamma(self.z[idx], self.x[idx], self.alpha_n, self.alpha_t) for idx in _indices)
        _print("γ =", _gamma)

        _d_amax = tuple(d_amax(self.z[idx], self.m_t, self.alpha_t, _gamma[idx]) for idx in _indices)
        _print("d_amax =", _d_amax)
        assert self.d_a[Ritzel] <= _d_amax[Ritzel]
        assert self.d_a[Rad] <= _d_amax[Rad]

        return