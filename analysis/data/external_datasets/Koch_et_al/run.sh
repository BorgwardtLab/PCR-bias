cd Bunny_P
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_M
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_F1
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_F2
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_F3
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_F4
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../Bunny_F5
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out