Contributing
============

Bug reports, feature suggestions and other contributions are greatly
appreciated! While I can't promise to implement everything, I will always try
to respond in a timely manner.

Short version
-------------

* Submit bug reports and feature requests at
  `GitHub <https://github.com/aburrell/apexpy/issues>`_
* Make pull requests to the ``develop`` branch

Bug reports
-----------

When `reporting a bug <https://github.com/aburrell/apexpy/issues>`_ please
include:

* Your operating system name and version
* Any details about your local setup that might be helpful in troubleshooting
* Detailed steps to reproduce the bug

Feature requests and feedback
-----------------------------

The best way to send feedback is to file an issue at
`GitHub <https://github.com/aburrell/apexpy/issues>`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions
  are welcome :)

Development
-----------

To set up `apexpy` for local development:

1. `Fork apexpy on GitHub <https://github.com/aburrell/apexpy/fork>`_.
2. Clone your fork locally::

    git clone git@github.com:your_name_here/apexpy.git

3. Create a branch for local development based off of the ``develop`` branch::

    git checkout -b name-of-your-bugfix-or-feature origin/develop

   Now you can make your changes locally. Add tests for bugs and new features
   in the relevant test file in the ``tests`` directory. The tests are run with
   ``pytest`` and can be written as normal functions (starting with ``test_``)
   containing a standard ``assert`` statement for testing output.

4. When you're done making changes, run ``pytest`` locally if you can::

    python -m pytest

5. Commit your changes and push your branch to GitHub::

    git add .
    git commit -m "ACRONYM: Brief description of your changes"
    git push origin name-of-your-bugfix-or-feature

   The project now uses the `NumPy acronyms <https://numpy.org/doc/stable/dev/development_workflow.html?highlight=development%20workflow>`_
   for development workflow in the commit messages.

6. Submit a pull request through the GitHub website. Pull requests should be
   made to the ``develop`` branch. The continuous integration (CI) testing
   servers will automatically test the whole codebase, including your changes,
   for multiple versions of Python on both Windows and Linux.

Pull Request Guidelines
^^^^^^^^^^^^^^^^^^^^^^^

If you need some code review or feedback while you're developing the code, just
make a pull request.

For merging, you should:

1. Include passing tests for your changes
2. Update/add documentation if relevant
3. Add a note to ``CHANGELOG.rst`` about the changes
4. Add yourself to ``AUTHORS.rst`` and ``.zenodo.json`` with your
   `ORCiD <https://orcid.org/>`_

Style Guidelines
^^^^^^^^^^^^^^^^

In general, apexpy follows PEP8 and numpydoc guidelines.  PyTest is used to run
the unit and integration tests, flake8 checks for style, and sphinx-build
performs documentation tests.  However, there are certain additional style
elements that have been settled on to ensure the project maintains a consistent
coding style:

- Line breaks should occur before a binary operator (ignoring flake8 W503)
- Preferably break long lines on open parentheses instead of using backslashes
- Use no more than 80 characters per line
- Several dependent packages have common nicknames, including:

  * ``import datetime as dt``
  * ``import numpy as np``

- Provide tests with informative failure statements and descriptive, one-line
  docstrings.

apexpy is working on modernizing its code style to adhere to these guidelines.
