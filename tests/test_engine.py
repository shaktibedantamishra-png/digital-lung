import sys, numpy as np, pandas as pd
sys.path.insert(0, 'engine')
import digital_lung as dl
P = dl.BASE


def test_1_mp_hand():
    # hand calculations (K = 1.5/30 = 0.05, R_aw = 10):
    # 0.098*20*(0.5*0.42^2/0.04 + 20*0.05*10*0.42^2 + 10*0.42) = 0.098*20*8.169 = 16.01124
    # 0.098*15*(0.5*0.5^2/0.05 + 15*0.05*10*0.5^2 + 5*0.5)     = 0.098*15*6.875 = 10.10625
    # 0.098*30*(0.5*0.3^2/0.03 + 30*0.05*10*0.3^2 + 15*0.3)    = 0.098*30*7.35  = 21.609
    for (C, vt, rr, pe), exp in [((0.04, 0.42, 20, 10), 16.01124), ((0.05, 0.5, 15, 5), 10.10625),
                                 ((0.03, 0.3, 30, 15), 21.609)]:
        assert abs(dl.mp(C, vt, rr, pe, P) - exp) < 1e-9


def test_2_matched_mv_closed_form():
    C, pbw, pe = 0.042, 70.0, 10.0
    K = (1 + P['IE']) / (60 * P['IE'])
    for v in dl.VT_GRID:
        rr = 100.0 / v; vt = v * pbw / 1000; VE = rr * vt
        closed = 0.098 * (0.5 * VE * vt / C + K * P['R_aw'] * VE**2 + pe * VE)
        assert abs(dl.mp(C, vt, rr, pe, P) - closed) < 1e-9


def test_3_fractions_sum():
    pe, pl = np.meshgrid(np.arange(0, 25.1, 0.5), np.arange(0, 60.1, 0.5))
    pl = np.maximum(pl, pe)
    fo, fc, fx = dl.fractions(pe, pl, P)
    assert np.allclose(fo + fc + fx, 1) and (fo >= 0).all() and (fc >= 0).all() and (fx >= -1e-12).all()


def test_4_pod_monotone_base():
    vts = np.arange(4, 10.01, 0.5); peeps = np.arange(5, 20.01, 0.5)
    for C in [0.015, 0.025, 0.04, 0.06, 0.09]:
        for pbw in [45.0, 60.0, 75.0, 90.0]:
            V, PE = np.meshgrid(vts * pbw / 1000, peeps, indexing='ij')
            _, od, _ = dl.penalties(C, V, 20.0, PE, P)
            assert (np.diff(od, axis=0) >= -1e-12).all(), ('vt', C, pbw)
            assert (np.diff(od, axis=1) >= -1e-12).all(), ('peep', C, pbw)


def test_5_pdr_proportional_rr():
    C, vt, pe = 0.03, 0.6, 5.0
    vals = [dl.penalties(C, vt, rr, pe, P)[0] / rr for rr in [10, 15, 20, 30, 45]]
    assert vals[0] > 0 and np.allclose(vals, vals[0])


def test_6_reproducible():
    pairs = pd.DataFrame({'db': np.repeat(dl.DBS, 50), 'crs': np.tile(np.linspace(20, 70, 50), 3),
                          'pbw': np.tile(np.linspace(45, 90, 50), 3)})
    a = dl.resample(pairs, n=100); b = dl.resample(pairs, n=100)
    assert a.equals(b)
    ca = dl.resample(pairs, n=100, method='copula'); cb = dl.resample(pairs, n=100, method='copula')
    assert ca.equals(cb)


if __name__ == '__main__':
    for k, f in list(globals().items()):
        if k.startswith('test_'): f(); print('PASS', k)
