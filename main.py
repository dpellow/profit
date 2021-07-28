#!/usr/bin/env python

import itertools
import numpy as np

import argparse

def parse_user_input():

    parser = argparse.ArgumentParser(
        description=
        'Output possible probe sets',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument('-t','--target',
     help='fasta file of the target sequences to find probes for',
     required=True, type=str
        )
    parser.add_argument('-n','--negative',
     help='fasta file of the negtive sequence set that probes cannot hit',
     required=True, type=str
        )
    parser.add_argument('-o','--outdir',
     help='output directory',
     required=True, type=str
        )
    parser.add_argument('-l','--probe_len',
     help='length of probes',
     required=False, type=int, default=20
        )
    parser.add_argument('-d','--max_degenerate',
     help='maximum number of degenerate bases allowed in a probe',
     required=False, type=int, default=1
        )
    parser.add_argument('-c','--min_coverage',
     help='minimum number of probes required to consider a target as covered',
     required=False, type=int, default=15
        )
    parser.add_argument('-f','--max_false',
     help='maximum number of times a negative sequence may be hit',
     required=False, type=int, default=3
        )
    return parser.parse_args()

def readfq(fp): # this is a generator function
    ''' Adapted from https://github.com/lh3/readfq
    '''
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def get_rc(seq):
    ''' Return the reverse complement of seq
    '''
    rev = reversed(seq)
    return "".join([complements.get(i,i) for i in rev])



def get_all_positions(infile,k, ndegen):
    """ Create a dictionary of all degenerate canonical kmers appearing in the sequences
        of infile mapped to the sequence and position of their appearance.
        Top level of dictionary is number of degenerate positions
    """

    pdict = {}

    masks = {i:[] for i in range(ndegen+1)}
    for i in masks:
        for c in itertools.combinations(range(k),i):
            masks[i].append(c)

    fp = open(infile)
    for name, seq, _ in readfq(fp): # Note: name of sequence can't have spaces in it
        l = len(seq)
        if l < k: continue
        for i in range(l-k+1):
            kmer = seq[i:i+k]
            rc = get_rc(kmer)
            if rc < kmer: kmer = rc
            if kmer not in pdict:
                pdict[kmer] = {'positive':{}, 'negative':{}}
            if name in pdict[kmer]['positive']:
                pdict[kmer]['positive'][name].append(i)
            else:
                pdict[kmer]['positive'][name] = [i]
    fp.close()
    return pdict


def get_negative_positions(pdict,infile,k):
    """ Fill in the positions of any kmer in pdict that appears in the negative set
    """
    fp = open(infile)
    for name, seq, _ in readfq(fp): # Note: name of sequence can't have spaces in it
        l = len(seq)
        if l < k: continue
        for i in range(l-k+1):
            kmer = seq[i:i+k]
            rc = get_rc(kmer)
            if rc < kmer: kmer = rc
            if kmer in pdict:
                if name in pdict[kmer]['negative']:
                    pdict[kmer]['negative'][name].append(i)
                else:
                    pdict[kmer]['negative'][name] = [i]
    fp.close()



def main(args):

    target_file = args.target
    neg_file = args.negative
    outdir = args.outdir
    plen = args.probe_len
    num_degen = args.max_degenerate
    num_probes = args.min_coverage
    num_fp = args.max_false

    MIN_DIST = 3 # minimum distance between probes - TODO: make this a parameter?

    # create a dictionary of positions of all probe length subsequences across all target sequences
    probes_dict = get_all_positions(target_file, plen, num_degen)

    # TODO: this assumes that the dictionary for the negative set will be too large,
    # i.e if it is all sequences from a huge metagenome.
    # if that assumption is not correct, then we can just build another dictionary for the negative
    # set instead of doing this.
    get_negative_positions(probes_dict,neg_file,plen)




if __name__=='__main__':
    args = parse_user_input()
    main(args)
