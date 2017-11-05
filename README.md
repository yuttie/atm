# atm
An implementation of the purity measure and more.


## Requirements
- C++ compiler
- Boost C++ Libraries
- CMake (to run the accompanying build script)

At least it should be able to be built with the following combinations:
- Gentoo Linux, {g++ 5.4, g++ 6.4, g++ 7.2, clang++ 5.0}, Boost 1.65, CMake 3.9.4
- Ubuntu 16.04, {g++ 5.4, clang++ 4.0}, Boost 1.58.0, CMake 3.5.1

The g++ compilers prior to the version 5.1 lacks necessary features and cannot be used.


## Build
Running the following commands yields a compiled program file under the directory `/path/to/atm/directory/build/Release`.

```sh
cd /path/to/atm/directory
cmake .
make
```

By default, this builds in Release mode, which results in an optimized program for performance with no debug information.
You can also specify a different "build type" when you run cmake as follows:

```sh
cmake . -DCMAKE_BUILD_TYPE=Debug
cmake . -DCMAKE_BUILD_TYPE=Release
cmake . -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake . -DCMAKE_BUILD_TYPE=MinSizeRel
```

According to the build type you specified, program files are placed under a directory `/path/to/atm/directory/build/<build-type>`.
If you don't know which build type to choose, you should go with the default Release build.


## Publications
- Yuta Taniguchi, Ryuji Masui, Toshihiro Aoyama, and Daisuke Ikeda,
  "Probabilistic Model for Purity Values of Bacterial Genome Sequences,"
  International Journal of Bioscience, Biochemistry and Bioinformatics, vol. 5, no. 5, pp. 288–295, Sep. 2015.
  https://doi.org/10.17706/ijbbb.2015.5.5.288-295
- Yuta Taniguchi, Yasuhiro Yamada, Osamu Maruyama, Satoru Kuhara, and Daisuke Ikeda,
  "The Purity Measure for Genomic Regions Leads to Horizontally Transferred Genes,"
  Journal of Bioinformatics and Computational Biology, vol. 11, no. 6, p. 1343002, Dec. 2013.
  https://doi.org/10.1142/S0219720013430026
