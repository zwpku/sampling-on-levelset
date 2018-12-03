# sampling-on-levelset

Illustrative codes for numerical sampling schemes on a level set.

These codes have been used in the paper 

   [Ergodic SDEs on submanifolds and related numerical sampling schemes](https://arxiv.org/abs/1702.08064)

## COMPILE & INSTALL

1. Download the source code.

```
	git clone https://github.com/zwpku/sampling-on-levelset.git
```

   	The code should be avaiable in the directory ./sampling-on-levelset

2. Enter the directory containing source files.

```
  	cd ./sampling-on-levelset/src
```

3. Compile.

```
    g++ gen_path.cpp com.c linpack.c ranlib.c -o gen_path
```

## DIRECTORIES & CODES
1. src (c/cpp codes)

     - gen_path.cpp:   implementation of different sampling schemes on an ellipse in R<sup>2</sup>.

     - com.c, linpack.c, ranlib.c, ranlib.h:    needed to generate random numbers, borrowed from the package [RANLIB.C](http://www.netlib.org/random/) (file: [ranlib.c.tar.gz](http://www.netlib.org/random/ranlib.c.tar.gz)).
     
2. plot-scripts (python scripts for generating figures)

	- elliptic.py:  illustration for the effects of different projection maps

	- plot_bin_from_traj.py: output density profiles for 3 different schemes

	- vec_field.py: plot the streamlines of different projection maps

	- path_check_scatter_plot.py: scatter plot for the states sampled using the Euler-Maruyama discretization

## USAGE
1.   Under the directory ./sampling-on-levelset, create directories for data and figures

```
    mkdir fig data
```

2.   Under the src directory, run the program gen_path, to generate data using different schemes. Data will be stored in ./sampling-on-levelset/data
     
3.   Run the scripts in  ./sampling-on-levelset/plot-scripts/ to plot figures. Figures will be generated in ./sampling-on-levelset/fig
