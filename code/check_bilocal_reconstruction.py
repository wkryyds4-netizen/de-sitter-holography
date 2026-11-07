"""Finite CAR and independent degenerate-bath tests of bilocal reconstruction.

Only this new main runs. The output is exclusive-create. These finite checks do
not establish a double-scaled limit or the asymptotic removal of the remainder.
"""
from datetime import datetime, timezone
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import platform
import time
import traceback

import numpy as np
from numpy.polynomial.legendre import leggauss
import scipy
from scipy.linalg import null_space

from complex_syk_dictionary import sample_hamiltonian, comm
from bilocal_channel import blocks

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / 'results/susskind-bilocal-reconstruction/finite-volterra-checks.json'
PARAMETERS = {
    'CAR': [{'N': 4, 'p': 2, 'm': 3, 'seed': 110701},
            {'N': 5, 'p': 2, 'm': 3, 'seed': 110702},
            {'N': 4, 'p': 2, 'm': 2, 'seed': 110703}],
    'generic': {'seed': 110704, 'h_eigenvalues': [-0.4, -0.4, 0.7],
                'bath_dimension': 4, 'alpha': 0.64, 'm': 3,
                'coupling_scale': 0.25},
    'time_tests': [0.37, 0.91], 'quadrature_orders': [8, 14],
    'contraction_times': [-1.4, -0.5, 0.0, 0.37, 0.91, 1.4],
    'laplace_z': [[0.7, 0.9], [1.3, -0.4]],
    'small_time_steps': [0.04, 0.02, 0.01],
    'pinching_eigenvalue_tolerance': 1e-10,
    'matrix_tolerance': 2e-11, 'quadrature_tolerance': 3e-10,
    'finite_nonzero_threshold': 1e-10,
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(a):
    a = np.asarray(a)
    return {'shape': list(a.shape), 'real': a.real.tolist(), 'imag': a.imag.tolist()}


def hs(a):
    return float(np.linalg.norm(a) / np.sqrt(a.shape[1]))


def rel(a, b):
    return float(np.linalg.norm(a-b) / max(1., np.linalg.norm(a), np.linalg.norm(b)))


class Compression:
    """Compressed spectral exponential, without computing a full exp each call."""
    def __init__(self, H, V):
        self.e, u = np.linalg.eigh(H)
        self.w = u.conj().T @ V

    @lru_cache(maxsize=None)
    def value(self, t):
        return self.w.conj().T @ (np.exp(1j*t*self.e)[:, None] * self.w)

    def derivative(self, t):
        return self.w.conj().T @ ((1j*self.e*np.exp(1j*t*self.e))[:, None]*self.w)

    def laplace(self, z):
        return self.w.conj().T @ ((1/(z-1j*self.e))[:, None]*self.w)


class Fixture:
    def __init__(self, h, a, bath, R, m, metadata, direct=None):
        self.energy, U = np.linalg.eigh(h)
        self.h = np.diag(self.energy)
        self.a = U.conj().T @ a @ U
        self.R = R @ U
        self.bath = bath
        self.m, self.d, self.metadata = m, len(h), metadata
        self.U, self.direct = U, direct
        self.groups = []
        for i, e in enumerate(self.energy):
            if not self.groups or abs(e-self.energy[self.groups[-1][0]]) > PARAMETERS['pinching_eigenvalue_tolerance']:
                self.groups.append([i])
            else:
                self.groups[-1].append(i)
        self.keep = np.zeros((self.d, self.d), dtype=bool)
        for ix in self.groups:
            self.keep[np.ix_(ix, ix)] = True
        self.memory = Compression(bath, self.R)
        K = np.block([[self.a, self.R.conj().T], [self.R, bath]])
        V = np.vstack([np.eye(self.d), np.zeros((len(bath), self.d))])
        self.full = Compression(K, V)
        self.bg_blocks = []
        for ix in self.groups:
            rg = self.R[:, ix]
            kg = np.block([[self.a[np.ix_(ix, ix)], rg.conj().T], [rg, bath]])
            vg = np.vstack([np.eye(len(ix)), np.zeros((len(bath), len(ix)))])
            self.bg_blocks.append((ix, Compression(kg, vg)))

    def pinch(self, x):
        return np.where(self.keep, x, 0.)

    def off(self, x):
        return np.where(self.keep, 0., x)

    @lru_cache(maxsize=None)
    def f0(self, t):
        out = np.zeros((self.d, self.d), complex)
        for ix, block in self.bg_blocks:
            out[np.ix_(ix, ix)] = block.value(t)
        return out

    def f0prime(self, t):
        out = np.zeros((self.d, self.d), complex)
        for ix, block in self.bg_blocks:
            out[np.ix_(ix, ix)] = block.derivative(t)
        return out

    @lru_cache(maxsize=None)
    def dm(self, t):
        return self.off(self.memory.value(t))

    def centered(self, x, t):
        phase = np.exp(-0.5j*t*self.energy)
        return phase[:, None]*x*phase[None, :]

    def field(self, t):
        return self.m*self.centered(self.dm(t), t)

    def bilocal(self, t):
        return self.m*self.off(self.centered(self.full.value(t), t))

    def matrices(self):
        return {k: encode(getattr(self, k)) for k in ['h', 'a', 'bath', 'R']}


def car_fixture(spec):
    n, p, m, seed = [spec[k] for k in ['N', 'p', 'm', 'seed']]
    cs, J, H = sample_hamiltonian(n, p, seed)
    scale = max(1., float(np.linalg.norm(H, 2)))
    H = H/scale
    inds, h, C, K = blocks(n, p, H, cs, m)
    V = C/np.sqrt(m)
    Q = null_space(V.conj().T)
    bath = Q.conj().T@K@Q
    R = Q.conj().T@K@V
    a = V.conj().T@K@V
    assert rel(a, (1-p/m)*h) < PARAMETERS['matrix_tolerance']
    # Full-Fock Heisenberg definition is independent of the reduced bath.
    ev, uv = np.linalg.eigh(H)
    def direct(t):
        u = (uv*np.exp(0.5j*t*ev))@uv.conj().T
        value = sum((u.conj().T@c@u).conj().T@(u@c@u.conj().T) for c in cs)
        return value[np.ix_(inds, inds)]
    fixture = Fixture(h, a, bath, R, m, {'kind': 'actual_complex_CAR', **spec,
                      'H_scale_divisor': scale, 'J_before_H_rescaling': encode(J),
                      'zero_memory_control': m == p}, direct)
    A = sum(comm(H, c).conj().T@comm(H, c) for c in cs)[np.ix_(inds, inds)]
    y = fixture.off(fixture.U.conj().T@A@fixture.U)
    assert rel(y, m*fixture.dm(0.)) < PARAMETERS['matrix_tolerance']
    fixture.metadata['force_initial_identity_error'] = rel(y, m*fixture.dm(0.))
    return fixture


def generic_fixture():
    spec = PARAMETERS['generic']
    rng = np.random.default_rng(spec['seed'])
    h = np.diag(spec['h_eigenvalues'])
    raw = rng.normal(size=(4, 4))+1j*rng.normal(size=(4, 4))
    bath = (raw+raw.conj().T)/2
    bath /= max(1., float(np.linalg.norm(bath, 2)))
    R = spec['coupling_scale']*(rng.normal(size=(4, 3))+1j*rng.normal(size=(4, 3)))
    return Fixture(h, spec['alpha']*h, bath, R, spec['m'],
                   {'kind': 'independent_degenerate_Hermitian_bath', **spec})


def gauss(order, length):
    x, w = leggauss(order)
    return (x+1)*length/2, w*length/2


def time_integrals(f, t, order):
    linear = np.zeros((f.d, f.d), complex)
    remainder = linear.copy()
    full_difference = linear.copy()
    filter_direct = linear.copy()
    filter_spectral = linear.copy()
    memory_ode = linear.copy()
    delta = f.energy[:, None]-f.energy[None, :]
    ss, ww = gauss(order, t)
    for s, ws in zip(ss, ww):
        dm = f.dm(float(s));L = t-s
        field = f.field(float(s))
        chi = L*np.sinc(delta*L/(2*np.pi))
        filter_spectral -= ws*chi*field
        memory_ode += ws*f.pinch(f.memory.value(float(t-s)))@f.f0(float(s))
        aa, wa = gauss(order, L)
        for a, w in zip(aa, wa):
            b = L-a
            left = f.f0(float(a))@dm
            f0b, fb = f.f0(float(b)), f.full.value(float(b))
            linear -= ws*w*(left@f0b)
            remainder -= ws*w*(left@(fb-f0b))
            full_difference -= ws*w*(left@fb)
            u = (a-b)/2
            filter_direct -= ws*w*np.exp(1j*u*delta)*field
    D = f.full.value(t)-f.f0(t)
    linear_b = f.m*f.off(f.centered(linear, t))
    remainder_b = f.m*f.off(f.centered(remainder, t))
    return {
        'order': order,
        'ordered_first_difference_error': rel(D, full_difference),
        'ordered_second_difference_error': rel(D, linear+remainder),
        'background_Volterra_error': rel(f.f0prime(t), 1j*f.a@f.f0(t)-memory_ode),
        'double_vs_spectral_filter_error': rel(filter_direct, filter_spectral),
        'bilocal_linear_plus_remainder_error': rel(f.bilocal(t), linear_b+remainder_b),
        'finite_remainder_HS': hs(remainder),
        'finite_projected_scaled_remainder_HS': hs(remainder_b),
        'finite_background_replacement_difference_HS': hs(linear_b-filter_spectral),
        'finite_actual_minus_free_filter_HS': hs(f.bilocal(t)-filter_spectral),
        '_linear': linear, '_remainder': remainder, '_filter': filter_spectral,
    }


def check_fixture(f):
    matrix_tol, quad_tol = PARAMETERS['matrix_tolerance'], PARAMETERS['quadrature_tolerance']
    contraction = []
    for t in PARAMETERS['contraction_times']:
        norm = float(np.linalg.norm(f.f0(t), 2))
        assert norm <= 1+matrix_tol
        assert hs(f.off(f.f0(t))) < matrix_tol
        contraction.append({'time': t, 'F0_operator_norm': norm,
                            'block_dilation_dimensions': [len(ix)+len(f.bath) for ix in f.groups]})
    resolvents = []
    for re, im in PARAMETERS['laplace_z']:
        z = complex(re, im)
        mhat = f.memory.laplace(z)
        m0 = f.pinch(mhat);dm = mhat-m0
        fhat = f.full.laplace(z)
        f0hat = np.zeros_like(fhat)
        for ix, block in f.bg_blocks:
            f0hat[np.ix_(ix, ix)] = block.laplace(z)
        expected0 = np.linalg.inv(z*np.eye(f.d)-1j*f.a+m0)
        exact_D = fhat-f0hat
        linear = -f0hat@dm@f0hat
        remainder = f0hat@dm@f0hat@dm@fhat
        checks = {'background_Schur_error': rel(f0hat, expected0),
                  'full_Schur_error': rel(fhat, np.linalg.inv(z*np.eye(f.d)-1j*f.a+mhat)),
                  'ordered_first_difference_error': rel(exact_D, -f0hat@dm@fhat),
                  'ordered_second_difference_error': rel(exact_D, linear+remainder)}
        assert max(checks.values()) < matrix_tol
        resolvents.append({'z': [re, im], **checks, 'finite_remainder_HS': hs(remainder),
                            'incorrect_order_difference_HS': hs(f0hat@dm@fhat-dm@f0hat@fhat)})
    times = []
    for t in PARAMETERS['time_tests']:
        field_identity = rel(f.dm(t), f.centered(f.field(t), -t)/f.m)
        assert field_identity < matrix_tol
        scans = [time_integrals(f, t, order) for order in PARAMETERS['quadrature_orders']]
        changes = {key: rel(scans[0][key], scans[1][key]) for key in ['_linear', '_remainder', '_filter']}
        for scan in scans:
            assert max(v for k, v in scan.items() if k.endswith('_error')) < quad_tol
            for key in ['_linear', '_remainder', '_filter']:
                del scan[key]
        assert max(changes.values()) < quad_tol
        adjoint_error = rel(f.bilocal(-t), f.bilocal(t).conj().T)
        background_removed = hs(f.off(f.centered(f.f0(t), t)))
        assert max(adjoint_error, background_removed) < matrix_tol
        row = {'time': t, 'memory_centering_identity_error': field_identity,
               'bilocal_adjoint_error': adjoint_error, 'background_after_pinching_HS': background_removed,
               'quadrature': scans, 'resolution_differences': changes}
        if f.direct is not None:
            direct = f.U.conj().T@f.direct(t)@f.U
            error = rel(direct, f.m*f.centered(f.full.value(t), t))
            assert error < matrix_tol
            row['full_Fock_bilocal_compression_error'] = error
        times.append(row)
    Y = f.m*f.dm(0.)
    small = []
    for dt in PARAMETERS['small_time_steps']:
        even = (f.bilocal(dt)+f.bilocal(-dt))/2
        curvature = -2*even/dt**2
        small.append({'time': dt, 'even_curvature_absolute_HS_error': hs(curvature-Y),
                      'even_curvature_relative_error': hs(curvature-Y)/max(hs(Y), 1e-14)})
    nonzero = hs(Y) > PARAMETERS['finite_nonzero_threshold']
    if nonzero:
        errors = [r['even_curvature_relative_error'] for r in small]
        assert errors[1] < errors[0]/3 and errors[2] < errors[1]/3
        assert max(r['finite_remainder_HS'] for r in resolvents) > PARAMETERS['finite_nonzero_threshold']
        assert max(r['incorrect_order_difference_HS'] for r in resolvents) > PARAMETERS['finite_nonzero_threshold']
    else:
        assert f.metadata.get('zero_memory_control')
        assert max(hs(f.memory.value(t)) for t in PARAMETERS['time_tests']) < matrix_tol
    block_commutator = hs(comm(f.pinch(f.memory.value(.3)), f.pinch(f.memory.value(.8))))
    diagonal_only_loss = hs(f.pinch(f.memory.value(.3))-np.diag(np.diag(f.memory.value(.3))))
    if f.metadata['kind'] == 'independent_degenerate_Hermitian_bath':
        assert [len(g) for g in f.groups] == [2, 1]
        assert block_commutator > 1e-5 and diagonal_only_loss > 1e-3
    return {'parameters': f.metadata, 'system_dimension': f.d, 'bath_dimension': len(f.bath),
            'energy_block_dimensions': [len(g) for g in f.groups], 'input_matrices': f.matrices(),
            'contraction_checks': contraction, 'resolvent_checks': resolvents, 'time_checks': times,
            'dynamic_initial_memory_HS': hs(Y), 'small_time_curvature': small,
            'degenerate_background_commutator_HS': block_commutator,
            'incorrect_diagonal_only_pinching_loss_HS': diagonal_only_loss}


def main():
    if OUT.exists():
        raise FileExistsError(OUT)
    started = time.perf_counter()
    records = []
    for fixture in [*(car_fixture(s) for s in PARAMETERS['CAR']), generic_fixture()]:
        record = check_fixture(fixture)
        records.append(record)
        print('checked', record['parameters']['kind'], record['system_dimension'],
              record['energy_block_dimensions'], flush=True)
    sources = ['code/check_bilocal_reconstruction.py', 'code/complex_syk_dictionary.py',
               'code/bilocal_channel.py', 'literature/susskind-memory-to-bilocal-sources.md']
    output = {'status': 'passed', 'created_utc': datetime.now(timezone.utc).isoformat(),
              'parameters': PARAMETERS, 'records': records,
              'counts': {'systems': len(records), 'actual_CAR_systems': 3, 'degenerate_bath_systems': 1,
                         'Laplace_points': 8, 'time_quadrature_runs': 16, 'contraction_times': 24,
                         'small_time_curvatures': 12},
              'elapsed_seconds': time.perf_counter()-started,
              'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                              'scipy': scipy.__version__},
              'input_sha256': {name: sha(ROOT/name) for name in sources},
              'scope': 'Finite deterministic matrices and floating-point quadrature with two resolutions. '
                       'No interval certification, no ensemble estimate and no DS extrapolation. '
                       'Nonzero finite remainders and background replacement errors are retained. '
                       'CAR uses standard anticommutators, normalized sector HS norms and the stated '
                       'rescaled string-time Hamiltonian. Pinching retains entire degenerate energy blocks.'}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open('x', encoding='utf-8') as stream:
        json.dump(output, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print('passed', output['counts'], 'elapsed', output['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        # Preserve the failure separately; never turn it into a successful record.
        directory = OUT.parent/'diagnostics'
        directory.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
        with (directory/f'finite-volterra-failure-{stamp}.json').open('x', encoding='utf-8') as stream:
            json.dump({'status': 'failed', 'parameters': PARAMETERS,
                       'script_sha256': sha(Path(__file__)), 'traceback': traceback.format_exc()}, stream, indent=2)
        raise
