PhD-python
==========

This repository has python modules I created to perform Fermi analysis during my PhD Work.

All of my python modules are under the 'lande' namespace. In order to use them, you need to clone this repository and add it to you rpython path:

```
$ git clone https://github.com/joshualande/PhD-python.git
$ cd PhD-python
$ export PYTHONPATH=$PWD
$ ipython -c 'from lande.fermi.likelihood import bandfitter'
```

All of the python modules can be found inside of the base lande package.
