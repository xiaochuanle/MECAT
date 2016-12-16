# <a name="S-installation"></a>Installation of MECAT-CA

* Step 1, download [wgs-8.2-Linux_amd64.tar.bz2](https://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.2/) and decompress it.
* Step2, copy `PBcR-Assembler` and `runCA` into `wgs-8.2/Linux-amd64/bin/` and replace the original ones.
* Step 3, go to folder `wgs-8.2/Linux-amd64/bin/` and run the following comands
```shell
chmod +x PBcR-Assembler & perl -e PBcR-Assembler
chmod +x runCA & perl -e PBcR-Assembler
```