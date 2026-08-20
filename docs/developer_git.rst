.. _Sec:developer-git:

Code architecture and Git workflow
==================================

This guide describes the structure of SORA and proposes a Git workflow for
changing the code, reviewing pull requests, and publishing releases. For
documentation dependencies, authoring, and build commands, see
:ref:`Sphinx documentation <Sec:developer-sphinx>`.

The examples assume that contributions target ``develop`` and that stable
releases are promoted to ``master``.

The terminal commands use `Git <https://git-scm.com/>`_. Forks, pull requests,
reviews, merges, and releases are managed in the GitHub web interface. Replace
usernames, branch names, pull request numbers, and versions with the values for
your work.


Create the code environment
---------------------------

SORA supports Python 3.11 or newer. Before running the Git commands, create and
activate an isolated Python environment with either ``venv`` or Conda.

Using ``venv`` (replace ``sora-dev`` if you prefer another directory name):

.. code-block:: console

   $ python -m venv sora-dev
   $ source sora-dev/bin/activate

On Windows PowerShell, activate the environment with
``sora-dev\Scripts\Activate.ps1``.

Using Conda (replace ``sora-dev`` if you prefer another environment name):

.. code-block:: console

   $ conda create --name sora-dev python=3.11 pip
   $ conda activate sora-dev

Keep the selected environment activated for the remaining steps.


Install Git
-----------

Git is required for the local workflow. Conda users can install it in the
active environment:

.. code-block:: console

   $ conda install --channel conda-forge git

If you are using ``venv``, install Git on the operating system. On Debian or
Ubuntu, use the system package manager:

.. code-block:: console

   $ sudo apt update
   $ sudo apt install git

On macOS with Homebrew:

.. code-block:: console

   $ brew install git

On Windows with WinGet, run these commands in PowerShell:

.. code-block:: powershell

   winget install --id Git.Git -e --source winget

Native installers may require a new terminal. If so, reactivate the selected
Python environment with the corresponding command from the previous section.
Then verify that Git is available:

.. code-block:: console

   $ git --version


Architecture at a glance
------------------------

SORA is an object-oriented Python package built around the physical components
of a stellar occultation. The package-level API in ``sora/__init__.py`` exposes
the main domain objects, while the subpackages contain their implementations
and supporting numerical functions.

::

    Star + Body
             \ Body owns an Ephem
              \
               + event time ----> prediction
                                  |
    Observer + LightCurve ------> Occultation
                                  |
                                  +--> Chord / ChordList
                                  |
                                  +--> fitting + stats
                                           |
                                           +--> size, shape and astrometry

The central objects cooperate as follows:

* :class:`sora.Star` represents the occulted star and handles catalogue data,
  astrometry, proper motion, and apparent diameter.
* :class:`sora.Body` represents the occulting Solar System body. It gathers
  physical data, shape and reference-frame information, and owns the ephemeris
  used to calculate its position.
* The ephemeris classes provide a common position interface backed by JPL
  Horizons, SPICE kernels, JPL files, or an ephemeris table.
* :class:`sora.Observer` and :class:`sora.Spacecraft` locate an observation in
  space, while :class:`sora.LightCurve` stores the photometry or already
  determined immersion and emersion times.
* :class:`sora.Occultation` is the integration point. It combines the star,
  body, ephemeris, and event time, then manages observations as
  :class:`sora.occultation.Chord` objects through a
  :class:`sora.occultation.ChordList`.
* ``sora.prediction`` searches for events and builds prediction tables and
  maps. ``sora.occultation.fitting`` and ``sora.stats`` turn projected chord
  positions into fitted physical results.


Repository layout
-----------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Path
     - Responsibility
   * - ``sora/__init__.py``
     - Public package facade; it re-exports the version resolved by
       ``sora/version.py``.
   * - ``sora/version.py``
     - Runtime version loader. A source checkout reads Git metadata through
       ``setuptools-scm``; installed distributions use the generated
       ``sora/_version.py`` module.
   * - ``sora/body/``
     - Solar System body data. ``core.py`` defines ``Body``;
       ``meta.py`` validates physical and ephemeris attributes; ``utils.py``
       contains database and magnitude helpers; ``frame/`` and ``shape/``
       implement body-fixed frames and projected shapes.
   * - ``sora/star/``
     - Stellar data, catalogue adapters, astrometric propagation, magnitudes,
       and angular diameter.
   * - ``sora/ephem/``
     - The ``BaseEphem`` contract and concrete Horizons, JPL, SPICE kernel, and
       table-backed ephemerides.
   * - ``sora/observer/``
     - Ground observers, spacecraft, and projected observer positions.
   * - ``sora/lightcurve/``
     - Light-curve input, normalization, occultation detection, and event
       modelling.
   * - ``sora/prediction/``
     - Event search, closest-approach parameters, prediction tables, and maps.
   * - ``sora/occultation/``
     - End-to-end event orchestration, chords, sky-plane positions, and limb or
       ellipse fitting.
   * - ``sora/stats/``
     - Parameters and numerical optimization backends used by fitting code.
   * - ``sora/extra/``
     - Shared plotting, ellipse, and chi-square utilities.
   * - ``sora/config/``
     - Cross-cutting input validation, deprecation helpers, collections, and
       terminal progress output.
   * - ``sora/data/``
     - Package data distributed with SORA.
   * - ``docs/``
     - Sphinx sources, API reference, examples, tutorials, and release notes.
       See :ref:`Sphinx documentation <Sec:developer-sphinx>`.
   * - ``pyproject.toml``
     - Build backend, package metadata, supported Python versions, runtime and
       optional dependencies, ``setuptools-scm``, pytest, and coverage
       configuration.
   * - ``setup.py``
     - Minimal compatibility shim that delegates the build to the declarative
       configuration in ``pyproject.toml``.
   * - ``MANIFEST.in``
     - Source-distribution inclusion and exclusion rules.
   * - ``tox.ini``
     - Isolated test environments for Python 3.11--3.13 and the documentation
       build environment.
   * - ``.github/workflows/ci.yml``
     - Distribution build and the OpenAstronomy reusable tox test matrix.


Following a call through the code
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the public facade first, then follow the object into its subpackage. For
example, ``from sora import Body`` resolves through ``sora/__init__.py`` and
``sora/body/__init__.py`` to ``sora/body/core.py``.

``Body`` is also a useful example of the internal design:

* ``Body.__init__`` validates keywords with ``sora.config.input_tests``.
* Data may come from the internal satellite data, JPL Small-Body Database, or
  explicit local values.
* ``BaseBody`` properties convert and validate physical values as
  ``PhysicalData`` quantities.
* Assigning an ephemeris shares the body identifiers and selected physical
  values with the ephemeris object.
* Frame and shape implementations are delegated to ``sora.body.frame`` and
  ``sora.body.shape``.

When changing an object, inspect its ``core.py``, base or metadata class,
utilities, package ``__init__.py``, API page under ``docs/modules/``, and
guideline notebook together. A public symbol should be exported deliberately;
an internal helper should remain private.


Configure a fork
----------------

Contributors without write access should create the fork in the browser. Follow
the `GitHub browser workflow for forks
<https://docs.github.com/en/pull-requests/how-tos/work-with-forks/fork-a-repo>`_:

#. Sign in to GitHub and open the `canonical SORA repository
   <https://github.com/riogroup/SORA>`_.
#. Click :guilabel:`Fork` in the upper-right corner.
#. Choose the owner for the fork and keep ``SORA`` as the repository name.
#. Choose whether to copy only the default branch. Additional branches can be
   fetched from ``upstream`` later.
#. Click :guilabel:`Create fork`.

If the fork already exists, open it in the browser and skip these steps. Clone
the fork locally in either case (replace ``YOUR-USERNAME`` with your GitHub
username):

.. code-block:: console

   $ git clone https://github.com/YOUR-USERNAME/SORA.git
   $ cd SORA

The clone configures your fork as ``origin``. Add the canonical repository as
``upstream``, then verify and fetch both remote references:

.. code-block:: console

   $ git remote add upstream https://github.com/riogroup/SORA.git
   $ git remote -v
   $ git fetch upstream --prune --tags

Both ``origin`` and ``upstream`` may contain a remote ``develop`` branch. On
the first checkout, name the canonical branch explicitly and create a local
tracking branch:

.. code-block:: console

   $ git switch --track upstream/develop

This creates the local ``develop`` branch with ``upstream/develop`` as its
upstream, avoiding ambiguity with ``origin/develop``. If local ``develop``
already exists, use ``git switch develop`` instead.

Maintainers who cloned the canonical repository as ``origin`` can use
``origin`` in place of ``upstream`` in the commands below.


Install SORA locally
--------------------

After cloning the fork and entering the repository root, ensure that the
environment created earlier is active, upgrade ``pip``, and install SORA in
editable mode with the test and documentation dependency groups declared in
``pyproject.toml``:

.. code-block:: console

   $ python -m pip install --upgrade pip
   $ python -m pip install -e ".[test,docs]"

Install the isolated-environment and distribution-build frontends as needed:

.. code-block:: console

   $ python -m pip install tox build

Some SORA operations query external services such as JPL Horizons, SBDB, MPC,
VizieR, and LIneA TAP. Keep network-dependent checks separate from deterministic
local checks, and record the service and input used when a change depends on an
online result.


Make a change
-------------

With the local tracking branch configured above, start every feature or fix
from an up-to-date ``develop``:

.. code-block:: console

   $ git fetch upstream --prune
   $ git switch develop
   $ git pull --ff-only upstream develop
   $ git switch -c feat/short-description

Suggested branch prefixes are ``feat/``, ``fix/``, ``docs/``, ``refactor/``,
and ``release/``. Keep one logical change per branch.

While working, inspect the exact patch and stage files intentionally:

.. code-block:: console

   $ git status --short
   $ git diff
   $ git add sora/path/to_changed_file.py docs/path/to_changed_page.rst
   $ git diff --cached
   $ git commit -m "Brief imperative summary"

Before opening the pull request, incorporate recent upstream changes. Merging is
the safest option for an already shared branch because it does not rewrite its
history:

.. code-block:: console

   $ git fetch upstream
   $ git merge upstream/develop
   $ git push -u origin feat/short-description


Implementation conventions
~~~~~~~~~~~~~~~~~~~~~~~~~~

Preserve the domain types used throughout the code: ``astropy.time.Time`` for
epochs, ``astropy.units.Quantity`` for dimensional values, and
``astropy.coordinates`` objects for positions and frames. State units,
coordinate direction, reference frame, and scalar-versus-vector behavior in
docstrings.

For a public API change:

#. update the implementation and its NumPy-style docstring;
#. review the subpackage ``__all__`` and ``__init__.py`` exports;
#. update the corresponding ``docs/modules/`` page or guideline;
#. add a release-note entry under the affected module;
#. preserve compatibility with a deprecation helper from
   ``sora.config.decorators`` when removal does not need to be immediate.

Keep catalogue and ephemeris service calls behind their existing adapters.
Numerical code should not acquire an unrelated network dependency.


Validate a change
-----------------

Run checks from the repository root. The baseline checks cover importability,
the source and wheel distributions, and the documentation:

.. code-block:: console

   $ python -m compileall -q sora
   $ python -c "import sora; print(sora.__version__)"
   $ python -m build
   $ tox run -e build_docs

Also run a focused example that exercises the changed behavior. Prefer local
inputs for the first check. When tracked tests are available, run the tox
environment for an installed interpreter (replace ``py311`` with ``py312`` or
``py313`` when appropriate):

.. code-block:: console

   $ tox run -e py311

The continuous-integration workflow always builds the source distribution and
wheel on Python 3.13. It enables the OpenAstronomy reusable tox workflow for
Python 3.11, 3.12, and 3.13 when Git contains a file matching
``sora/**/tests/test_*.py``. A local untracked test does not enable that matrix;
add intentional tests to the pull request.

When a change affects docstrings, tutorials, release notes, or the public API,
also follow :ref:`the Sphinx validation workflow <Sec:developer-sphinx>`.

For code changes, verify at least:

* expected output and failure behavior;
* units, time scales, and coordinate frames;
* scalar and array inputs when the API accepts both;
* offline or service-error behavior for network-backed features;
* backward compatibility for an existing public interface.


.. _Sec:developer-git-open-pr:

Open a pull request
-------------------

Push the branch to your fork:

.. code-block:: console

   $ git push -u origin feat/short-description

Then create the pull request in the browser, following GitHub's `pull request
from a fork workflow
<https://docs.github.com/en/pull-requests/how-tos/create-pull-requests/creating-a-pull-request-from-a-fork>`_:

#. Open the `canonical SORA repository <https://github.com/riogroup/SORA>`_
   and select :guilabel:`Pull requests`.
#. Click :guilabel:`New pull request`, then :guilabel:`compare across forks`.
#. Select ``riogroup/SORA`` and ``develop`` as the base repository and branch.
#. Select your fork as the head fork and ``feat/short-description`` as the
   compare branch.
#. Review the commits and changed files, then enter a concise title and a full
   description.
#. Optionally enable :guilabel:`Allow edits from maintainers` after considering
   the additional warning GitHub displays for forks containing Actions
   workflows.
#. Click :guilabel:`Create pull request`, or use the adjacent menu to create a
   draft when the change is not ready for review.

The pull request description should explain the problem, the chosen solution,
scientific or API implications, validation commands and results, external
services used, and documentation or release-note changes. Keep unrelated
cleanup out of the pull request.

After receiving feedback, add normal commits to the same branch and push again.
The pull request updates automatically:

.. code-block:: console

   $ git add path/to/file
   $ git commit -m "Address review feedback"
   $ git push


Review a pull request
---------------------

Begin in the browser with the stated goal, metadata, patch, and checks:

#. Open the canonical repository and click :guilabel:`Pull requests`.
#. Select the pull request to review.
#. Use :guilabel:`Conversation` to read its goal, discussion, linked issues,
   reviewers, and merge status.
#. Inspect :guilabel:`Commits` and :guilabel:`Files changed`, then check the
   reported workflow results before running the change locally.

For local validation, first make sure your current work is committed or safely
stored. Find the pull request number in the browser, then fetch its GitHub
reference into a temporary local branch (replace ``123`` as needed):

.. code-block:: console

   $ git status --short
   $ git fetch upstream pull/123/head:review/pr-123
   $ git switch review/pr-123
   $ python -m build

If the pull request includes tracked tests, run the applicable tox environment
as well:

.. code-block:: console

   $ tox run -e py311

Apply the focused checks described by the author. A scientific review should
also trace units and frames, test representative edge cases, and distinguish
deterministic calculations from remote-service behavior. Review dependency and
public API changes explicitly. If the pull request changes documentation, build
it as described in :ref:`Sphinx documentation <Sec:developer-sphinx>`.

Submit the review through GitHub's `browser review workflow
<https://docs.github.com/en/pull-requests/how-tos/review-pull-requests/reviewing-proposed-changes-in-a-pull-request>`_:

#. Open :guilabel:`Files changed` and use the line comment button for focused
   feedback or suggestions. Start a review when comments should be grouped.
#. Click :guilabel:`Review changes` and write a concise summary.
#. Select exactly one outcome: :guilabel:`Comment` for non-blocking feedback,
   :guilabel:`Approve` when the change is ready, or
   :guilabel:`Request changes` for a blocking issue.
#. Click :guilabel:`Submit review`.

Authors cannot approve their own pull requests. Maintainers should merge only
after required reviews and checks are satisfied. In the pull request's
:guilabel:`Conversation` tab, scroll to the merge box, select the repository's
intended merge method, click its merge button, and confirm. Optionally delete
the source branch afterward. See GitHub's `browser merge procedure
<https://docs.github.com/en/pull-requests/how-tos/merge-and-close-pull-requests/merging-a-pull-request>`_.


Publish a release
-----------------

Release operations require maintainer access to the canonical GitHub
repository and the configured package index accounts. Use a release branch to
make the promotion reviewable.

Prepare the release
~~~~~~~~~~~~~~~~~~~

Create the branch from the current ``develop``:

.. code-block:: console

   $ git fetch upstream --prune --tags
   $ git switch develop
   $ git pull --ff-only upstream develop
   $ git switch -c release/v0.3.4

Create ``docs/releases/v0.3.4.rst`` with the release notes and date, then add it
at the top of ``docs/releases.rst``. Summarize new features, API changes, bug
fixes, and documentation, and link the relevant pull requests or issues with
the existing Sphinx roles.

Do not edit a version string in ``setup.py``, ``sora/__init__.py``, or
``docs/conf.py``. The project declares a dynamic version in ``pyproject.toml``;
``setuptools-scm`` derives it from Git tags, ``sora.__version__`` exposes it at
runtime, and Sphinx imports that value. The generated ``sora/_version.py`` file
is ignored in a source checkout and must not be committed.

Build and inspect the candidate distribution:

.. code-block:: console

   $ python -m pip install build twine
   $ python -m compileall -q sora
   $ python -m build --outdir build/release-candidate
   $ python -m twine check build/release-candidate/*

Before the release tag exists, ``setuptools-scm`` gives these candidate files a
development version derived from the previous tag and current commit. Check the
exact ``0.3.4`` filenames only after creating the release tag.

Before committing the release, complete the build described in
:ref:`Sphinx documentation <Sec:developer-sphinx>`.

Commit the release notes and push the branch:

.. code-block:: console

   $ git add docs/releases.rst docs/releases/v0.3.4.rst
   $ git commit -m "Prepare release v0.3.4"
   $ git push -u origin release/v0.3.4

Open the pull request in the browser with the process described in
:ref:`Open a pull request <Sec:developer-git-open-pr>`. Select ``master`` as
the base branch and ``release/v0.3.4`` from your fork as the compare branch.
Use a title such as ``Release v0.3.4`` and explain that the pull request
promotes the validated development line and prepares the release.

Tag and create the GitHub release
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After the release pull request is approved and merged, create an annotated tag
on the exact merge commit on ``master``. The local tag makes
``setuptools-scm`` resolve the final version. Build and validate the artifacts
with that tag present before pushing it:

.. code-block:: console

   $ git fetch upstream --prune
   $ git switch master
   $ git pull --ff-only upstream master
   $ git tag -a v0.3.4 -m "SORA v0.3.4"
   $ python -c "import sora; print(sora.__version__)"
   $ python -m build --outdir dist/v0.3.4
   $ python -m twine check dist/v0.3.4/sora_astro-0.3.4*
   $ git push upstream v0.3.4

Create the GitHub release in the browser, following the `GitHub release
procedure
<https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository>`_:

#. Open the canonical repository and click :guilabel:`Releases`.
#. Click :guilabel:`Draft a new release`.
#. Open :guilabel:`Choose a tag` and select the existing ``v0.3.4`` tag. Do not
   create a replacement tag in the form.
#. Confirm that the selected tag is the one pushed from the reviewed
   ``master`` commit.
#. Select the previous release tag and click :guilabel:`Generate release notes`,
   or write the notes manually from the reviewed release-note file.
#. Set the title to ``SORA v0.3.4`` and, if desired, attach the already checked
   files from ``dist/v0.3.4``.
#. Review the draft, then click :guilabel:`Publish release`. Use
   :guilabel:`Save draft` instead if another maintainer must verify it first.

Publish the Python package
~~~~~~~~~~~~~~~~~~~~~~~~~~

Test the previously checked artifacts on TestPyPI before publishing the same
files to PyPI:

.. code-block:: console

   $ python -m twine upload --repository testpypi \
       dist/v0.3.4/sora_astro-0.3.4*
   $ python -m pip install --index-url https://test.pypi.org/simple/ \
       --extra-index-url https://pypi.org/simple/ sora-astro==0.3.4
   $ python -m twine upload dist/v0.3.4/sora_astro-0.3.4*

Never rebuild artifacts between TestPyPI and PyPI: testing and publishing the
same files ensures that the verified distribution is the released
distribution. Credentials must come from the maintainer's configured keyring,
token, or trusted publishing environment and must never be committed.

Synchronize ``master`` back into ``develop`` after the release so the release
notes and tagged commit remain in the development line. Make the sync through a
pull request:

.. code-block:: console

   $ git switch -c sync/v0.3.4-into-develop upstream/master
   $ git push -u origin sync/v0.3.4-into-develop

Open another pull request in the browser. Select ``develop`` as the base branch
and ``sync/v0.3.4-into-develop`` from your fork as the compare branch, then use
``Sync v0.3.4 into develop`` as the title. Review and merge it through the same
browser workflow used for other pull requests.


Release checklist
~~~~~~~~~~~~~~~~~

* ``develop`` contains the intended features and fixes.
* The release tag is the single version source, and its value agrees with the
  runtime API, documentation, and built artifacts.
* Release notes describe user-visible changes and link their pull requests.
* Tox tests (when present), documentation, import smoke check, focused behavior
  checks, and distribution validation pass.
* The release pull request targets ``master`` and is approved.
* The annotated tag points to the reviewed commit on ``master``.
* The GitHub release and package index use the same version.
* The exact TestPyPI artifacts are uploaded to PyPI.
* ``master`` is synchronized back into ``develop``.
