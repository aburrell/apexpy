#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import sys


if __name__ == "__main__":
    str_to_bool = {"false": False, "true": True}

    base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    print("Project path: {0}".format(base_path))

    env_path = os.path.join(base_path, ".tox", "bootstrap")
    if sys.platform == "win32":
        bin_path = os.path.join(env_path, "Scripts")
    else:
        bin_path = os.path.join(env_path, "bin")

    if not os.path.exists(env_path):
        print("Making bootstrap env in: {0} ...".format(env_path))
        try:
            subprocess.check_call(["virtualenv", env_path], shell=True)
        except Exception:
            subprocess.check_call([sys.executable, "-m", "virtualenv",
                                   env_path], shell=True)

        print("Installing `jinja2` and `matrix` into bootstrap environment ...")
        subprocess.check_call([os.path.join(bin_path, "pip"), "install",
                               "jinja2", "matrix"], shell=True)

    activate = os.path.join(bin_path, "activate_this.py")
    act_dict = {__file__: activate}

    # tox requires activation with exec
    exec(compile(open(activate, "rb").read(), activate, "exec"), act_dict)

    # Import here to ensure success, since they may have just been installed
    import jinja2
    import matrix

    jinja = jinja2.Environment(
        loader=jinja2.FileSystemLoader(os.path.join(base_path, "ci",
                                                    "templates")),
        trim_blocks=True, lstrip_blocks=True, keep_trailing_newline=True)
    tox_environments = {}
    setup_matrix = matrix.from_file(os.path.join(base_path,
                                                 "setup.cfg")).items()
    for (alias, conf) in setup_matrix:
        python = conf["python_versions"]
        deps = conf["dependencies"] if "dependencies" in conf.keys() else ""
        if "coverage_flags" in conf.keys():
            cover = str_to_bool[conf["coverage_flags"].lower()]
        if "environment_variables" in conf.keys():
            env_vars = conf["environment_variables"]

        tox_environments[alias] = {
            "python": "python" + python if "py" not in python else python,
            "deps": deps.split()}

        if "coverage_flags" in conf.keys():
            tox_environments[alias].update(cover=cover)
        if "environment_variables" in conf.keys():
            tox_environments[alias].update(env_vars=env_vars.split())

    for name in os.listdir(os.path.join("ci", "templates")):
        with open(os.path.join(base_path, name), "w") as fh:
            fh.write(jinja.get_template(name).render(
                tox_environments=tox_environments))
        print("Wrote {}".format(name))
    print("DONE.")
