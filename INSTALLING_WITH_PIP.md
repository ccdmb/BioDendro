# Installing python packages

BioDendro is a python program which can be used as a python module or via a command-line interface.
It can be installed from PyPI <https://pypi.org/project/BioDendro/> using [pip](https://pip.pypa.io/en/stable/).

It is tested to work with python 3.6+, and it depends on several packages:

- [numpy](http://www.numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/)
- [xlrd](https://github.com/python-excel/xlrd)
- [xlsxwriter](https://xlsxwriter.readthedocs.io/)
- [plotly](https://plot.ly/)


The install procedures below will install these for you automatically.


## Typical install with pip

```bash
python3 -m pip install --user BioDendro
```

The `--user` flag tells pip to install to a user directory rather than a system directory (e.g. somewhere under `/usr/` or `Program Files`).
To install as root or in Windows you can omit the `--user` flag, but you may need root or administrator permission.

Note that the `--user` subdirectory containing scripts (`<userdir>/bin` on linux/mac or `<userdir>\Scripts` on Windows), may not automatically be on your path.
The user directory installed to is given by the python command `import site; print(site.USER_BASE)`.
Generally on linux this is `~/.local`.
You can change this directory by setting the `$PYTHONUSERBASE` environment variable.
You can then add the `Script` or `bin` subdirectory to your `$PATH` environment variable
(e.g. [Linux and MacOS](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix), or [Windows](https://stackoverflow.com/questions/9546324/adding-directory-to-path-environment-variable-in-windows), or just google it ;) ).

Power users and sys-admins may also be interested in the `pip install --prefix` flag, but this will also require you to update your `$PYTHONPATH` environment variable or use [`.pth` files](https://docs.python.org/3/library/site.html) .
Nonetheless, it's useful for [module](http://modules.sourceforge.net/) software management.


## Installing in a `virtualenv`

Generally it's a good idea to install python packages in a python [virtual environment](https://virtualenv.pypa.io/en/stable/)
This protects versions of packages required for your operating system to work, from being updated (potentially breaking older code).

Here's a basic rundown of the `virtualenv` workflow.

```bash
# If it isn't installed already run one of these
# Try to use the system package managers if possible to avoid mixing up system dependencies.
# Mac users may be interested in homebrew.
sudo python3 -m pip install virtualenv
sudo apt install python3-virtualenv # Ubuntu and probably Debian
sudo dnf install python3-virtualenv # Fedora 24+

# Change dir to where you want the envronment to live (usually a project dir).
cd my_project

# Create a virtualenv in a folder ./env
# python3 can be substituted with your version of python.
python3 -m venv env

# Loads the virtualenv (essentially this changes PYTHONPATH and some other environment variables).
source env/bin/activate

which python3
# ./env/bin/python3

python3 -m pip install BioDendro

BioDendro --help
```

This time we don't need to use sudo or tell pip to install as `--user` because it will install all files in the `./env` folder.
Note however, that whenever you start a new terminal and want to use this environment you'll need to repeat the `source env/bin/activate`.
If you want to exit the virtual environment, you can type `deactivate` on the command line or just start a new terminal.

## Installing from git repository.

To install most recent commit on the `master` branch:

```bash
python3 -m venv env
source env/bin/activate
python3 -m pip install git+https://github.com/ccdmb/BioDendro.git
```

Note that we installed the package in a virtual environment, which is always recommended.
See the section above for more details.

You can also clone the repository yourself and then install from that.

```bash
git clone https://github.com/ccdmb/BioDendro.git ./BioDendro
cd BioDendro

python3 -m venv env
source env/bin/activate

python3 -m pip install .
```

This is particularly useful if you want to modify the code.
You can install the module in a way that allows changes you make to be automatically installed.

```bash
# Use this after cloning.
python3 -m venv env
source env/bin/activate

python3 -m pip install -e ".[dev,test]"
```

The `[dev,test]` thing just installs some extra packages used for testing etc.
