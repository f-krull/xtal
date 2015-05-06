# xtal
Collection of C++ tools to process PDB files. See src/xtal* subfolders for specific notes. 

compile: 
 ```
git clone https://github.com/f-krull/xtal
cd xtal
make
 ```

compile whithout gcc-ar:
 ```
make GCC_AR=ar
 ```

compile on Raspberry Pi:
 ```
make raspi GCC_AR=ar
 ```

## dependencies
* g++ >= 4.5 (for -flto)
* make
* wget
