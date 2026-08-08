# First try _dev.scm_version when setuptools-scm is available. The _dev
# package is excluded from wheels and tarballs, which then use the generated
# _version module instead.
try:
    try:
        from ._dev.scm_version import version
    except ImportError:
        from ._version import version
except Exception:
    import warnings

    warnings.warn(
        f'could not determine {__name__.split(".")[0]} package version; this indicates a broken installation'
    )
    del warnings

    version = '0.0.0'
