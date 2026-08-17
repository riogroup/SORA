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
   * - ``docs/requirements.txt``
     - Versioned Python dependencies needed to build the documentation.
   * - ``docs/modules/``
     - API reference generated from package docstrings with autodoc.
   * - ``docs/guidelines/``
     - Task-oriented Jupyter notebook tutorials and their input/output assets.
   * - ``docs/examples/``
     - Additional notebooks and example pages.
   * - ``docs/releases.rst``
     - Release-note index; version pages live under ``docs/releases/``.
   * - ``.readthedocs.yaml``
     - Read the Docs operating system, Python version, build command, and
       output formats.


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

Install SORA in editable mode and then install the versioned documentation
toolchain:

.. code-block:: console

   $ python -m pip install -e .
   $ python -m pip install -r docs/requirements.txt

The editable SORA installation is required because autodoc imports the package
and its runtime dependencies while building the API reference.
``docs/requirements.txt`` installs Sphinx, ``nbsphinx``, the Read the Docs
theme, and the extensions used for issues, bibliography, and videos.


Install Pandoc
~~~~~~~~~~~~~~

The documentation contains Jupyter notebooks, so a full build also requires the
external `Pandoc <https://pandoc.org/installing.html>`_ executable on ``PATH``.
Pandoc is not a Python dependency and therefore is not installed by
``docs/requirements.txt``.

Conda users can prepare the complete environment with:

.. code-block:: console

   $ conda activate sora-dev
   $ python -m pip install -e .
   $ python -m pip install -r docs/requirements.txt
   $ conda install -c conda-forge pandoc

Do not use ``pip install pandoc`` for this requirement: the PyPI project with
that name does not provide the command-line executable that ``nbsphinx`` calls.

Verify that Python, Sphinx, and Pandoc are available to the same process:

.. code-block:: console

   $ command -v python
   $ python -m sphinx --version
   $ command -v pandoc
   $ pandoc --version


Build the HTML documentation
----------------------------

Run Sphinx from the repository root:

.. code-block:: console

   $ python -m sphinx -b html docs docs/build/html

Open ``docs/build/html/index.html`` in a browser and inspect every page changed
by the contribution. The ``docs/build`` directory is generated, excluded by
Sphinx, and ignored by Git.

Use a fresh Sphinx environment when cached source information may be stale:

.. code-block:: console

   $ python -m sphinx -E -b html docs docs/build/html

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

   $ conda run -n sora-dev python -m sphinx -b html docs docs/build/html

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

* SORA is installed in editable mode.
* Packages from ``docs/requirements.txt`` are installed.
* The real Pandoc executable is available on ``PATH``.
* The page is included in a ``toctree``.
* API changes update their docstrings and module reference.
* Notebook cells and asset paths are reproducible.
* The complete HTML build succeeds.
* Changed pages render correctly and introduce no new warnings.
