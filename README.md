## Please see the [wiki](https://github.com/hmdlab/RaptRanker/wiki) for detail usage

# RaptRanker

RaptRanker is a software for RNA aptamer selection from HT-SELEX experiment data based on local sequence and structure information.

If you have any issues or recommendations, please feel free to open a [ticket](https://github.com/hmdlab/HTAptamerSelection_VIP/issues).

## Requirements
RaptRanker needs C++ compiler, CMake, and Boost Libraries. 
We checked under the following versions.


- OS : CentOS 7.5.1804, macOS Mojave 10.14.5
- C++ compiler : g++ 5.5.0, clang 10.0.1
- CMake : 3.13.4, 3.14.2
- Boost : 1.69.0, 1.71.0

## Installation
Download or clone the latest version sources, and build with CMake.

```
git clone https://github.com/hmdlab/RaptRanker.git
cd RaptRanker
mkdir build
cd build
cmake ..
make
```

`$ warning: This header is deprecated. Use <boost/integer/integer_log2.hpp> instead.` may occur during make. This warning occurs in third-party source code. For now RaptRanker works, so please ignore it.

## Usage
The executable file will be complied as `RaptRanker/bin/RaptRanker`. RaptRanker input parameters with a json file. To work RaptRanker, create a parameter file and call
```
path/to/RaptRanker/bin/RaptRanker parameter.json
```

For the detail of parameters, please check [wiki](https://github.com/hmdlab/RaptRanker/wiki/The-parameter-file).
