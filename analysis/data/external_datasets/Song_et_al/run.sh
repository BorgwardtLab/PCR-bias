cd PCR1
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../PCR2
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../PCR3
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../PCR4
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../PCR5
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out

cd ../PCR6
~/.local/bin/bbmap/bbduk.sh -Xmx6g in1=R1.fq.gz in2=R2.fq.gz out=stdout.fq ref=adapters ktrim=r k=23 mink=11 hdist=1 tpe tbo 2> log_filtering.out | ~/.local/bin/bbmap/bbmap.sh -Xmx6g in=stdin.fq ref=design_files.fasta ordered interleaved nodisk outu=unmapped.fq.gz scafstats=scafstats.txt 2> log_mapping.out