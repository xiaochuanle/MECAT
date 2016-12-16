#</a>Contents
* [Introduction](#S-introduction)
* [Installation](#S-installation)
* [Quick Start](#S-quick-start)
* [Program Descriptions](#S-program-description)
* [Citation](#S-citation)
* [Contact](#S-contact)

# <a name="S-introduction"></a>Introdction

MECAT is an ultra-fast Mapping, Error Correction and de novo Assembly Tools for single molecula sequencing (SMRT) reads. MECAT employs novel alignment and error correction algorithms that are much more efficient than the state of art of aligners and error correction tools. MECAT can be used for effectively de novo assemblying large genomes. For example, on a 32-thread computer, MECAT takes 14 days to assemble a human genome based on 54x SMRT data, which is 40 times faster than the current [PBcR-Mhap pipeline](http://cbcb.umd.edu/software/pbcr/mhap/). We also use MECAT to assemble a diploid human genome based on 108x SMRT data only in 40 days. The latter assembly leads a great improvement of quality to the previous genome assembled from the 54x haploid SMRT data. MECAT performance were compared with [PBcR-Mhap pipeline](http://cbcb.umd.edu/software/pbcr/mhap/) and [FALCON](https://github.com/PacificBiosciences/falcon) in five real datasets. The quality of assembled contigs produced by MECAT is the same or better than that of the [PBcR-Mhap pipeline](http://cbcb.umd.edu/software/pbcr/mhap/) and [FALCON](https://github.com/PacificBiosciences/falcon). Here are some comparisons:

<div>
    <table border="0">
      <tr>
        <th>one</th>
        <th>two</th>
      </tr>
      <tr>
        <td>Hello</td>
        <td>你好</td>
      </tr>
    </table>
</div>

MECAT consists of four modules:
* `mecat2pw`, a fast and accurate pairwise mapping tool for SMRT reads
* `mecat2ref`, a fast and accurate reference mapping tool for SMRT reads
* `mecat2cns`, correct noisy reads based on their pairwise overlaps
* `mecat2canu`, a modified and more efficient version of the [Canu pipeline](https://github.com/marbl/canu). [Canu](https://github.com/marbl/canu) is a customized version of the Celera Assembler that designed for high-noise single-molecule sequencing

MEAP is written in C, C++, and perl. It is open source and distributed under the [GPLv3](http://www.gnu.org/licenses/gpl-3.0.html) license.

# <a name="S-installation"></a>Installation

```
git clone https://github.com/marbl/canu.git
cd mecat
./intall.sh
```

After installation, all the executables are found in `mecat/Linux-amd64/bin`. The folder name `Linux-amd64` will vary in operating systems. For example, in MAC, the executables are put in `mecat/Darwin-amd64/bin`.

# <a name="S-quick-start"></a>Quick Start

Using MECAT to assemble a genome involves 4 steps. Here we take assemblying the genome of Ecoli as an example, to go through each step in order. Options of each command will be given in next section. We download the reads [ecoli_filtered.fastq.gz](http://gembox.cbcb.umd.edu/mhap/raw/ecoli_filtered.fastq.gz) from the MHAP website. By decompressing it we obtain `ecoli_filtered.fastq`.

Before assemblying, we must add `mecat/Linux-amd64/bin` to the system path. For example,
```shell
export PATH=/public/users/chenying/mecat/Linux-amd64/bin:$PATH
```
And we assume that all the following commands take place in the same folder.

* Step 1, using `mecat2pw` to detect overlapping candidates
``` shell
mecat2pw -j 0 -d ecoli_filtered.fastq -o ecoli_filtered.fastq.pm.can -w wrk_dir -t 16
```

* Step 2, correct the noisy reads based on their pairwise overlapping candidates.
```shell
mecat2cns -i 0 -t 16 ecoli_filtered.fastq.pm.can ecoli_filtered.fastq corrected_ecoli_filtered.fasta
```

* Step 3, extract the longest 25X corrected reads
```shell
extract_sequences corrected_ecoli_filtered.fasta corrected_ecoli_25x.fasta 4800000 25
```

* Step 4, assemble the longest 25X corrected reads using `mecat2cacu`
```shell
mecat2canu -assemble -p ecoli -d ecoli genomeSize=4800000 ErrorRate=0.02 maxMemory=40 maxThreads=16 useGrid=0 Overlapper=mecat2asmpw -pacbio-corrected corrected_ecoli_25x.fasta
```

After step 4, the assembled genome is given in file `ecoli/ecoli.contigs.fasta`. Details of the contigs can be found in file `ecoli/ecoli.layout.tigInfo`.

# <a name="S-program-description"></a>Program Descriptions

We describe in detail each module of MECAT, including their options and output formats.

## <a name="SS-pairwise-mapping"></a>`mecat2pw`

### </a>options
The command for running `mecat2pw` is
```shell
mecat2pw -j [task] -d [fasta/fastq] -w [working folder] -t [# of threads] -o [output] -n [# of candidates] -g [0/1]
```
The options are:
* `-j [task]`, job name, 0 = detect overlapping candidates only, 1 = output overlaps in M4 format, default = 1. If we are to correct noisy reads, outputing overlapping candidates is enough.
* `-d [fasta/fastq]`, reads file name in FASTA or FASTQ format.
* `-w [working folder]`, a directory for storing temporary results, will be created if not exists.
* `-t [# of threads]`, number of CPU threads used for overlapping, default=1.
* `-o [output]`, output file name
* `-n [# of candidates]`, number of candidates considered for gapped extension, default=100. Since each chunk  is about 2GB size, number of candidates(NC) should be set by genome size (GS).For GS < 20M, NC should be set as 200; For GS>20M and GS<200M; NC should be set as 100; For GS>200M, NC should be set as 50.
* `-g [0/1]`, output the gapped extension start point (1) or not (0), default=0.

### </a>output format

If the job is detecting overlapping candidates, the results are output in `can` format, each result of which occupies one line and 9 fields:
```shell
[A ID] [B ID] [A strand] [B strand] [A gapped start] [B gapped start] [voting score] [A length] [B length]
```

`mecat2pw` outputs overlapping results in `M4` format, of which one result is given in one line. The fileds of `M4` format is given in order below:
```shell
[A ID] [B ID] [% identity] [voting score] [A strand] [A start] [A end] [A length] [B strand] [B start] [B end] [B length]
```
If the `-g` option is set to `1`, two more fields indicating the extension starting points are given:
```shell
[A ID] [B ID] [% identity] [voting score] [A strand] [A start] [A end] [A length] [B strand] [B start] [B end] [B length] [A ext start] [B ext start]
```
In the `strand` field, `0` stands for the forward strand and `1` stands for the reverse strand. All the positions are zero-based and are based on the forward strand, whatever which strand the sequence is mapped. Here are some examples:
```shell
44 500 83.6617 30 0 349 8616 24525 0 1 10081 21813
353 500 83.2585 28 0 10273 18410 22390 1 0 10025 21813
271 500 80.4192 13 0 14308 19585 22770 1 4547 10281 21813
327 501 89.8652 117 0 10002 22529 22529 1 9403 21810 21811
328 501 90.8777 93 0 0 10945 22521 1 0 10902 21811
```
In the examples above, read 500 overlaps with reads 44, 353, 271, 327 and 328.

### </a> memory consumption

Before overlapping is conducted, the reads will be split into several chunks. Each chunk is about 2GB in size so that the overlapping can be run on a 8GB RAM computer.

## <a name="SS-reference-mapping"></a> `mecat2ref`

###</a>options

`mecat2ref` is used for mapping SMRT reads to the reference genomes. The command  is
```shell
mecat2ref -d [reads] -r [reference] -w [folder] -t [# of threads] -o [output] -b [# of results] -m [output format]
```
The meanings of each option are as follows:
* `-d [reads]`, reads file name in FASTA/FASTQ format
* `-r [reference]`, reference genome file name in FASTA format
* `-w [folder]`, a directory for storing temporary results
* `-t [# of threads]`, number of working CPU threads
* `-o [output]`, output file name
* `-b [# of result]`, output the best b alignments
* `-m [output format]`, output format: 0 = ref, 1 = M4, 2 = SAM, default = 0

###</a>output format

`mecat2ref` outputs results in one of the three formats: the `ref` format, the `M4` format, and the `SAM` format.

For the `ref` format, each result occupies three lines in the form:
```shell
[read name] [ref name] [ref strand] [voting score] [read start] [read end] [read length] [ref start] [ref end]
mapped read subsequence
mapped reference subsequence
```
The strands of the reads are always forward. In the `[ref strand]` field, `F` indicates forward strand while `R` indicates reverse strand. All the positions are zero-based and relative to the forward strand. Here is an example:
```shell
1	gi|556503834|ref|NC_000913.3|	F	10	2	58	1988134	1988197
AAT-AGCGCCTGCCAGGCG-TCTTTT--CCGGCCATTGT-CGCAG--CACTGTAACGCGTAAAA
AATTAGCGCCTGCCAGGCGGTCTTTTTTCCGGCCATTGTTCGCAGGG-ACTGTAACGCGTAAAA
```
In this example, read 1 is mapped to the reference `gi|556503834|ref|NC_000913.3|`.

###</a>memory consumption

* Index for the genome: genomeSize * 8 bytes
* Compressed index for each CPU thread: genomeSize \* 0.1 * t bytes
* Local alignment: 100M * t + 1G bytes

##<a name="SS-consensus"></a> `mecat2cns`

`mecat2cns` is an adaptive error correction tool for high-noise single-molecula sequencing reads. It is as accurate as [pbdagcon](https://github.com/PacificBiosciences/pbdagcon.git) and as fast as FalconSense. Inputs to `mecat2cns` can be either `can` format or `M4` format. The command for running `mecat2cns` is
```shell
mecat2cns [options] overlaps-file reads output
```
The options are
* `-i [input type]`, input format, 0 = `can`, 1 = `M4`, default = 1
* `-t [# of threads]`, number of CPU threads for consensus, default=1
* `-p [batch size]`, batch size the reads will be partitioned, default=100000
* `-r [ratio]`, minimum mapping ratio, default=0.9
* `-c [coverage]`, minimum coverage, default=6
* `-l [length]`, minimum length of the corrected sequence

If the inputs are `M4` format, the overlap results in `[overlaps-file]` must contain the gapped extension start point, which means the option `-g` in `mecat2pw` must be set to 1, otherwise `mecat2cns` will fail to run. Also note that the memory requirement of `mecat2cns` is about 1/4 of the total size of the reads. For example, if the reads are of total size 1GB, then `mecat2cns` will occupy about 250MB memory.

### </a> output format

The corrected sequences are given in FASTA format. The header of each corrected sequence consists of three components seperated by underlines:
```shell
>A_B_C_D
```
where
* `A` is the original read id
* `B` is the left-most effective position
* `C` is the right-most effective position
* `D` is the length of the corrected sequence

by effective position we mean the position in the original sequence that is covered by at least `c` (the argument to the option `-c`) reads.


##<a name="SS-assembly"></a>`mecat2canu`

`mecat2canu` is a modified and more efficient version of the [Canu pipeline](https://github.com/marbl/canu). `mecat2canu` accelerates `canu` by replacing its overlapper `mhap` by `mecat2asmpw`, which is a customized version of `mecat2pw`. The options of `mecat2canu` are the same as `canu` except its overlapper is replaced by `mecat2asmpw`. The command is
``` shell
mecat2canu -d [working-folder] -p [file-prefix] -assemble errorRate=[fraction error] \
	-overlapper=mecat2asmpw genomeSize=[genome size] \
    maxMemory=[host memory size] maxThreads=[# of CPU threads] usedGrid=0 \
    -pacbio-corrected reads-name
```

After assembling, the results are given in the folder `working-folder`. The assembled genome is given in the file `working-folder/file-prefix.contigs.fasta` and the details of the contigs are given in the file `working-folder/file-prefix.layout.tigInfo`.

#<a name="citation"></a>Citation
Chuan-Le Xiao, Ying Chen, Shang-Qian Xie, Kai-Ning Chen, Yan Wang, Feng Luo, Zhi Xie. [MECAT: an ultra-fast mapping, error correction and de novo assembly tool for single-molecule sequencing reads](http://biorxiv.org/content/early/2016/12/15/089250). bioRxiv 089250; doi: https://doi.org/10.1101/089250
#<a name="contact"></a> Contact
* Chuan-Le Xiao, xiaochuanle@126.com
* Ying Chen, chenying2016@gmail.com
* Feng Luo, luofeng@clemson.edu









