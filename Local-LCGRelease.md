## Installation of LCG dev3 in a local machine

- LCGCMake requires CMake installed in your system. On CVMFS (CERN) set `PATH` to use one of latest CMake versions (even 3.6.X is too old, LCGCMake does require CMake > 3.8). This is the only step that requires access to CVMFS:

```
ARCH=$(uname -m)
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.8.1/Linux-${ARCH}/bin:${PATH}
```

- Clone the LCGCMake repo:

```
mkdir $PWD/lcgcmake
git clone https://gitlab.cern.ch/sft/lcgcmake.git lcgcmake
```

- For convenience, let's create an alias to the main script:

```
alias lcgcmake=$PWD/lcgcmake/bin/lcgcmake
```

- Create a new directory to install ROOT and all their dependencies:

```
mkdir $PWD/localdev3
```

- Bootstrap the LCGCMake configuration:

```
lcgcmake configure --compiler=gcc62 --prefix=$PWD/localdev3 --version=dev3
```

- Build and install ROOT with its depedencies:

```
lcgcmake install ROOT
```

- After `ROOT` should be installed. To set the running environment do:

```
source $PWD/localdev3/dev3/x86_64-centos7-gcc62-opt/setup.sh
```

**NOTE** that in this last step your path may be different after `dev3` if you are using a different platform.

- Finally, we can check our installation running RDataFrame from Python, it should look like this:

```
$ python
Python 2.7.15 (default, Jul 30 2018, 09:33:44)
[GCC 6.2.0] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import ROOT
>>> df = ROOT.ROOT.RDataFrame(1)
```
