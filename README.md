# RamAnalyser

RamAnalyser is a simple tool to analyse the overlaps found by ram(https://github.com/lbcb-sci/ram)

# Usage

To build RamAnalyser run the following commands:

```bash
git clone https://github.com/Andrewzhang217/RamAnalyser && cd RamAnalyser && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. && make
```

which will create the RamAnalyser executable. Running the executable will display the following usage:

```bash
usage: ./RamAnalyser [options ...] <sequences>

  <target>/<sequences>
    input file in FASTA/FASTQ format (can be compressed with gzip)

  options:
    -k <std::uint8_t>
      default: 15
      length of minimizers
    -w <std::uint8_t>
      default: 5
      length of sliding window from which minimizers are sampled
    -f <double>
      default: 0.001
      threshold for ignoring most frequent minimizers
    -M <bool>
      default: false
      set minhash to true, use only a portion of all minimizers
    -t <std::uint8_t>
      default: 10
      number of threads
    -s <int>
      default: 1000
      sample size of sequences
    -h
      prints the usage
```
