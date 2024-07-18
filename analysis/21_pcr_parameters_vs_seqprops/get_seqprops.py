import argparse
import sys
import Bio.SeqIO
import subprocess
import os
import tempfile
import re
import tqdm.auto

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def get_longest_hp(seq):
    lchar = ''
    longest = 0
    for i in seq:
        if lchar == i:
            cnt += 1
        else:
            cnt = 1
        if cnt > longest:
            longest = cnt
        lchar = i
    return longest

pattern = re.compile('Minimum folding energy is (.*) kcal/mol.')
def get_folding_dg(seq):
    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        with open('tmp.fa', 'w') as f:
            f.write(f">tmp\n{seq}\n")
        p = subprocess.run(['mfold', "SEQ=tmp.fa", "NA=DNA"], capture_output=True)
        output = p.stdout.decode()
        match = pattern.search(output)
        if match:
            try:
                # print(f"Found minimum folding energy: {match.group(1)}")
                return float(match.group(1))
            except ValueError:
                print(f"Could not convert {match.group(1)} to float.")
                return None
        else:
            print(f"Could not find minimum folding energy in output: {output}")
            return None




def fasta2seqproperties(args):
    parser = argparse.ArgumentParser(description='Extracting sequence properties from a fasta file.')
    parser.add_argument('input', type=str, help='Path to the input file')
    parser.add_argument('output', type=str, help='Path to the output file')
    config = parser.parse_args(args)

    get_props = lambda seq: [
        seq.id,
        len(seq),
        seq.count('A')/len(seq),
        seq.count('C')/len(seq),
        seq.count('G')/len(seq),
        seq.count('T')/len(seq),
        (seq.count('G') + seq.count('C'))/len(seq),
        str(seq[0]),
        str(seq[-1]),
        get_longest_hp(str(seq.seq)),
        get_folding_dg(str(seq.seq)),
    ]
    
    i = 0
    with open(config.output, 'w') as fo:
        fo.write("id,length,A,C,G,T,GC,first,last,hp,dg\n")
        for seq in tqdm.auto.tqdm(Bio.SeqIO.parse(config.input, 'fasta')): 
            i += 1
            fo.write(f"{','.join(map(str, get_props(seq)))}\n")
    logger.info(f"Analysed {i} sequences from {config.input} and saved properties to {config.output}.")


if __name__ == "__main__":
    fasta2seqproperties(sys.argv[1:])