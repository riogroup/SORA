import numpy as np
import pytest

from sora.extra import ChiSquare


@pytest.mark.parametrize('invalid', [np.nan, np.inf, -np.inf])
def test_get_nsigma_excludes_nonfinite_trials(invalid):
    result = ChiSquare(
        np.array([1.0, invalid, 0.5, 3.0]),
        npts=4,
        immersion=np.array([10.0, 999.0, 10.1, 10.9]),
    )
    original = result.data['chi2'].copy()

    interval = result.get_nsigma()

    assert interval['chi2_min'] == 0.5
    assert interval['n_points'] == 2
    assert interval['sigma'] == 1
    np.testing.assert_allclose(interval['immersion'], [10.05, 0.05])
    np.testing.assert_array_equal(result.data['chi2'], original)
    np.testing.assert_allclose(result.get_nsigma(key='immersion'), [10.05, 0.05])


@pytest.mark.parametrize('chi2', [
    [np.nan, np.nan],
    [np.inf, -np.inf],
    [np.nan, np.inf, -np.inf],
    [],
])
def test_get_nsigma_requires_finite_trials(chi2):
    result = ChiSquare(np.array(chi2), npts=4)

    with pytest.raises(ValueError, match='No finite chi-square values'):
        result.get_nsigma()


@pytest.mark.parametrize('sigma', [1, 2, 3])
def test_get_nsigma_preserves_finite_confidence_selection(sigma):
    # Values exactly on the threshold remain excluded by the strict inequality.
    chi2 = np.array([0.5, 1.0, 1.5, 3.0, 4.5, 9.5])
    immersion = np.array([10.1, 10.0, 10.2, 10.9, 11.0, 12.0])
    result = ChiSquare(chi2, npts=8, immersion=immersion)
    selected = immersion[chi2 < chi2.min() + sigma ** 2]

    interval = result.get_nsigma(sigma=sigma)

    assert interval['chi2_min'] == chi2.min()
    assert interval['sigma'] == sigma
    assert interval['n_points'] == len(selected)
    np.testing.assert_allclose(interval['immersion'], [
        (selected.max() + selected.min()) / 2,
        (selected.max() - selected.min()) / 2,
    ])


def test_get_nsigma_single_finite_trial():
    result = ChiSquare(
        np.array([np.nan, 0.5, -np.inf]),
        npts=4,
        immersion=np.array([999.0, 10.1, -999.0]),
    )

    interval = result.get_nsigma()

    assert interval['n_points'] == 1
    assert interval['immersion'] == [10.1, 0.0]


def test_get_nsigma_rejects_unknown_parameter():
    result = ChiSquare(
        np.array([np.nan, 0.5]), npts=4, immersion=np.array([9.0, 10.1]),
    )

    with pytest.raises(ValueError, match='not one of the available keys'):
        result.get_nsigma(key='unknown')
