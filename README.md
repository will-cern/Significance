# Significance
Implementation of Formulae for estimating significance of event counts wrt a prediction with or without uncertainty

## Setting up the code

The recommended way to use this code is with one of the `StatAnalysis` releases where this package comes precompiled. Just setup like this:

```
asetup StatAnalysis,<version>
```

where `<version>` is e.g. 0.0.2.

### compile from source
Alternatively you can compile from source if you have a ROOT build available. You can either take the header and source files and compile it yourself, or use the provided `CMakeLists.txt` to compile and install like so:

```
mkdir build; cd build
cmake ../
make install
```

A library will be built in the `build` directory, which you can incorporate into your projects. By installing the library with the above target you can also use the code from the ROOT prompt.

## Using the code

The code can be used from the ROOT prompt or in ROOT macros with:

```
double n = 4; double predicted = 3; double uncert = 0.2;
Significance::Recommended(n, predicted, uncert);
```

Similarly, the code can be used from python with `PyROOT`:

```
import ROOT
ROOT.Significance.Recommended(4,3,0.2)
```


## Other formulae for significance

The code includes implementations of all the significance formulae that were investigated. The recommended method is accessed by the `Recommended` method, other methods can be accessed for curious developers by looking at the header.


