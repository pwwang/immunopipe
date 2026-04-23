import os

__version__ = "2.4.4"

desc = (
    f"Immunopipe v{__version__}\n"
    "Integrative analysis for scRNA-seq and scTCR-/scBCR-seq data"
)

_host_version = os.environ.get("IMMUNOPIPE_HOST_VERSION", None)
if _host_version and _host_version != __version__:
    desc = (
        f"{desc}\n\n"
        "WARNING: Immunopipe version mismatch:\n"
        f"  host version is 'v{_host_version}', but VM version is 'v{__version__}'.\n"
        "  This may cause unexpected errors. \n"
        "  Make sure to use the same version of immunopipe on host and VM."
    )
