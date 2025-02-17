# WisePIFinder

We propose a novel sketch algorithm called WisePIFinder, which aims to detect PI flows more accurately and efficiently in realtime. 

# How to run

Suppose you've already cloned the repository.

You just need:

```
$ make
$ ./WisePIFinder (-d dataset -m memory -k k)
```

**optional** arguments:

- -d: set the path of dataset to run, default dataset is CAIDA "1.dat"
- -m: set the memory size (KB), default memory is 100KB

# Output format

Our program will print the Throughput of insertion, AAE, ARE and Precision of these algorithms on the screen.
