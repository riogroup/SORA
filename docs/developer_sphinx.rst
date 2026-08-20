.. _Sec:developer-sphinx:

Sphinx documentation
====================

This guide explains how to install the documentation dependencies, edit the
Sphinx sources, build the site locally, and diagnose common build problems.
For the SORA architecture, branches, pull requests, and releases, see
:ref:`Code architecture and Git workflow <Sec:developer-git>`.


Documentation structure
-----------------------

The documentation is built from the ``docs/`` directory and published with
Read the Docs.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Path
     - Responsibility
   * - ``docs/index.rst``
     - Root document and main navigation trees.
   * - ``docs/conf.py``
     - Sphinx extensions, theme, project metadata, external references, and
       excluded build paths.
   * - ``pyproject.toml``
     - The ``docs`` optional-dependency group containing Sphinx, ``nbsphinx``,
       the theme, and the extensions used by the documentation.
   * - ``tox.ini``
     - The reproducible ``build_docs`` environment and its Sphinx command.
   * - ``docs/Makefile`` and ``docs/make.bat``
     - Convenience entry points for local Sphinx builds on Unix-like systems
       and Windows.
   * - ``docs/modules/``
     - API reference generated from package docstrings with autodoc.
   * - ``docs/guidelines/``
     - Task-oriented Jupyter notebook tutorials and their input/output assets.
   * - ``docs/examples/``
     - Additional notebooks and example pages.
   * - ``docs/releases.rst``
     - Release-note index; version pages live under ``docs/releases/``.
   * - ``.readthedocs.yaml``
     - Read the Docs operating system, Python version, ``docs`` extra,
       system-level Pandoc package, builder, and output formats.


Install the build dependencies
------------------------------

SORA supports Python 3.11 or newer. From the repository root, create and
activate an isolated environment with either ``venv`` or Conda.

Using ``venv``:

.. code-block:: console

   $ python -m venv .venv
   $ source .venv/bin/activate
   $ python -m pip install --upgrade pip

On Windows PowerShell, activate it with
``.venv\Scripts\Activate.ps1``.

Using Conda (replace ``sora-dev`` if you prefer another environment name):

.. code-block:: console

   $ conda create --name sora-dev python=3.11 pip
   $ conda activate sora-dev
   $ python -m pip install --upgrade pip

Install SORA in editable mode with the documentation dependencies declared by
the ``docs`` optional-dependency group in ``pyproject.toml``:

.. code-block:: console

   $ python -m pip install -e ".[docs]"

For direct Sphinx builds, the editable SORA installation is required because
autodoc imports the package and its runtime dependencies while building the API
reference. The ``docs`` extra installs Sphinx, ``nbsphinx``, the Read the Docs
theme, and the extensions used for issues, bibliography, and videos.

The canonical build uses tox, which creates an isolated environment and
installs the same ``docs`` extra automatically. Install tox in the active
environment with:

.. code-block:: console

   $ python -m pip install tox


Install Pandoc
~~~~~~~~~~~~~~

The documentation contains Jupyter notebooks, so a full build also requires the
external `Pandoc <https://pandoc.org/installing.html>`_ executable on ``PATH``.
Pandoc is not a Python dependency and therefore is not installed by the
``docs`` extra.

Conda users can prepare the complete environment with:

.. code-block:: console

   $ conda activate sora-dev
   $ python -m pip install -e ".[docs]"
   $ python -m pip install tox
   $ conda install -c conda-forge pandoc

Do not use ``pip install pandoc`` for this requirement: the PyPI project with
that name does not provide the command-line executable that ``nbsphinx`` calls.

Read the Docs handles the two dependency types separately according to
``.readthedocs.yaml``: it installs SORA with the ``docs`` extra and installs the
Pandoc executable from the operating system's APT packages. Local builds must
provide Pandoc with Conda or the appropriate system package manager.

Verify that Python, Sphinx, and Pandoc are available to the same process:

.. code-block:: console

   $ command -v python
   $ python -m sphinx --version
   $ command -v pandoc
   $ pandoc --version


Build the HTML documentation
----------------------------

Run the reproducible tox environment from the repository root:

.. code-block:: console

   $ tox run -e build_docs

Tox installs the project with the ``docs`` extra and invokes Sphinx with the
settings in ``tox.ini``. Open ``docs/_build/html/index.html`` in a browser and
inspect every page changed by the contribution. The ``docs/_build`` directory
is generated, excluded by Sphinx, and ignored by Git.

For a faster build in an already prepared editable environment, invoke Sphinx
directly:

.. code-block:: console

   $ python -m sphinx -b html docs docs/_build/html

The platform-specific wrappers are equivalent alternatives when their commands
are available:

.. code-block:: console

   $ make -C docs html

On Windows Command Prompt or PowerShell, run:

.. code-block:: powershell

   docs\make.bat html

Use a fresh Sphinx environment when cached source information may be stale:

.. code-block:: console

   $ python -m sphinx -E -b html docs docs/_build/html

The build should complete without introducing new warnings. Some existing
notebooks may emit warnings about old Pygments lexers or undefined internal
labels; compare the result with the ``develop`` baseline and investigate every
new warning related to the changed files.


Resolve ``PandocMissing``
~~~~~~~~~~~~~~~~~~~~~~~~~

Selecting a Conda environment's Python executable in an IDE does not always
activate that environment's ``PATH``. In that situation, Sphinx imports from
the expected environment while ``nbsphinx`` cannot find the adjacent Pandoc
executable.

Confirm the problem with:

.. code-block:: console

   $ python -c "import shutil; print(shutil.which('pandoc'))"

If this prints ``None``, activate the environment before building or run Sphinx
through Conda:

.. code-block:: console

   $ conda run -n sora-dev python -m sphinx -b html docs docs/_build/html

This command ensures that the Python modules and external executables come from
the same environment.


Edit Sphinx sources
-------------------

Use reStructuredText files for regular pages. Match the existing heading
hierarchy and leave blank lines around directives, lists, and code blocks.
Add every new page to a ``toctree`` so Sphinx can include it in navigation.

Use cross-references instead of hard-coded generated HTML paths:

.. code-block:: rst

   .. _Sec:example:

   Example title
   =============

   See :ref:`the example <Sec:example>`.

Public Python APIs are documented from their NumPy-style docstrings. When
changing a public class, function, parameter, return value, unit, or exception,
update both the docstring and the appropriate page under ``docs/modules/``.

Release notes are included from ``docs/releases.rst``. Add a new version file
under ``docs/releases/`` and include it at the top of the release index.
Use the configured ``:pr:``, ``:issue:``, and bibliography roles rather than
embedding repository URLs repeatedly.


Work with notebooks
~~~~~~~~~~~~~~~~~~~

``nbsphinx`` converts the notebooks under ``docs/guidelines/`` and
``docs/examples/``. Before committing a changed notebook:

* restart its kernel and execute cells in order when execution is safe;
* remove accidental tracebacks and machine-specific paths;
* keep downloaded or generated files out of Git unless they are intentional
  documentation assets;
* verify that referenced images and input files use paths relative to the
  notebook;
* rebuild the complete Sphinx site, not only the notebook itself.

Examples that depend on remote astronomical services should state the service
and input used. A temporary service failure should be distinguishable from an
error in the documented SORA behavior.


Documentation checklist
-----------------------

* Tox is installed for the canonical ``build_docs`` environment, or SORA is
  installed in editable mode with the ``docs`` extra for a direct build.
* The real Pandoc executable is available on ``PATH``.
* The page is included in a ``toctree``.
* API changes update their docstrings and module reference.
* Notebook cells and asset paths are reproducible.
* The complete HTML build succeeds.
* Changed pages render correctly and introduce no new warnings.
