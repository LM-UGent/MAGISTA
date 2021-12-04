#!/usr/bin/env python

# Copyright 2020 Department of Computational Biology for Infection Research - Helmholtz Centre for Infection Research
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import argparse
import os
#import gzip
import sys

class contig:
    id=""
    seq=""
    q=""

def fastx_as_dict(f):
    res= {}
    curentry = None
    mode= 0
    curid = 0
    for line in f:
        line=line.strip()
        if mode!=1 and (line[0]=='>' or line[0]=='@'):
             if curentry is not None :
                  res[curid] = curentry
                  curid += 1
             curentry = contig()
             curentry.id = line[1:]
             curentry.seq=""
             curentry.q = ""
             mode = 0
        elif mode==0:
             if line[0]=='+':
                 mode += 1
             else:
                 curentry.seq += line
        else:
             curentry.q += line
             if len(curentry.q) == len(curentry):
                 mode=0
    if curentry is not None:
        res[curid] = curentry
        curid += 1
    return [ res[x] for x in range(curid) ]



def read_lengths_from_fastx_file(fastx_file):
    """

    @param fastx_file: file path
    @type fastx_file: str
    @rtype: dict[str, int]
    """
    if ".gz" in fastx_file:
        raise RuntimeError("cannot handle gzip") 
    else:
        f = open(fastx_file, "rt")

    length = {}
    if os.path.getsize(fastx_file) == 0:
        return length

    for seq_record in fastx_as_dict(f):
        length[seq_record.id] = len(seq_record.seq)

    f.close()

    return length


def add_column(mapping_file, fasta_file):
    lengths = read_lengths_from_fastx_file(fasta_file)
    sequence_id_column_index = 0
    with open(mapping_file, 'r') as read_handler:
        for line in read_handler:
            line = line.rstrip('\n')
            if line.startswith('@@'):
                for index, column_name in enumerate(line[2:].split('\t')):
                    if column_name == 'SEQUENCEID':
                        sequence_id_column_index = index
                print('%s\t_LENGTH' % line)
                continue

            if len(line.strip()) == 0 or line.startswith('#') or line.startswith('@'):
                print(line)
                continue

            row_data = line.split('\t')
            sequence_id = row_data[sequence_id_column_index]
            if(sequence_id in lengths):
                print('%s\t%d' % (line, lengths[sequence_id]))


def main():
    parser = argparse.ArgumentParser(description="Add length column _LENGTH to gold standard mapping and print mapping on the standard output")
    parser.add_argument("-g", "--gold_standard_file", help="GSF", required=True)
    parser.add_argument("-f", "--fasta_file", help="FASTA or FASTQ file with sequences of gold standard", required=True)
    args = parser.parse_args()
    add_column(args.gold_standard_file, args.fasta_file)


if __name__ == "__main__":
    main()
