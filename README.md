# RDataFrame-Totem
---

## How to run `distill.py`

**Requirements**

This script reads Totem data from `eos`, namely from the following path:

```
/eos/totem/data/cmstotem/2015/90m/Totem/Ntuple/version2/4495/
```

 Therefore, the totem project needs to be mounted and accessible for the user.

### Using pure python from a terminal

1. Clone this repository:

```
git clone https://github.com/JavierCVilla/RDataFrame-Totem.git
```

2. Prepare the environment

  - The code requires `ROOT-6.14.00` or greater and `Python`.
  - Simplest way to fulfil this software dependencies is using the LCG Releases available through CVMFS.
  - The following command will setup your environment with these packages ready to be used:

  ```
  source /cvmfs/sft.cern.ch/lcg/views/dev3python3/latest/x86_64-slc6-gcc62-opt/setup.sh
  ```

  - Alternatively, your own `ROOT` and `Python` installation can be used, in which case you should ensure the python `ROOT` module is properly configured in your environment so it can be imported:

  ```
  python
  >>> import ROOT
  >>>
  ```

  - If the previous import failed, your `PYTHONPATH` may not be properly set. The easiest way to configure the environment for `root` is using its own setup script:

  ```
  source /your/path/to/root/bin/thisroot.sh
  ```

3. Run the code:

```
python distill.py <diagonal> [threads number]
```  
Valid diagonals: d45b_56t, d45t_56b, ad45b_56b, ad45t_56t

### Using HelixNebula

- Init a session in Swan HelixNebula and select the bleeding edge software stack, this is the only one that currently provides `ROOT-6.14.00`.
- Just copy `distill.ipynb` to your Cernbox space or to your SWAN instance in HelixNebula.
- Open the python notebook and execute the cells.
- HelixNebula already provides the needed environment configuration as well as access to the `eos` files.


## Comparing results

ROOT files produced by this code and the original analysis can be compared using the `rootcompare.c` script to ensure the same results are produced.

After setting up the environment, compile the script using:

```
g++ -o rootcompare rootcompare.c `root-config --cflags --glibs`
```

This program receives 4 arguments:

```
./rootcompare fileA treenameA fileB treenameB
```

**NOTE: This script is not meant to be generic enough to compare any pair of root files, currently it's aim to compare only files produced by this analysis.**
