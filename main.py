#!/usr/bin/env python

import bisect
import itertools
import numpy as np
from heapq import merge

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




def get_all_positions(infile,k, masks):
    """ Create a dictionary of all degenerate canonical kmers appearing in the sequences
        of infile mapped to the sequence and position of their appearance.
        Top level of dictionary is number of degenerate positions
    """

    ###{i:[] for i in range(ndegen+1)}
####    for i in masks:
        ####for c in itertools.combinations(range(k),i):
            ####masks[i].append(c)

    pdict = {m:{} for m in masks}
    targets = set()

    fp = open(infile)
    for name, seq, _ in readfq(fp): # Note: name of sequence can't have spaces in it
        l = len(seq)
        if l < k: continue
        targets.add(name)
        for i in range(l-k+1):
            kmer = seq[i:i+k]
            rc = get_rc(kmer)
            if rc < kmer: kmer = rc

            # all possible positions for degenerate nucleotides
            for m in masks:
                mkmer = ''.join([kmer[ind] for ind in range(k) if ind not in m])

                if mkmer not in pdict[m]:
                    pdict[m][mkmer] = {'positive':{}, 'negative':{}}
                if name in pdict[m][mkmer]['positive']:
                    pdict[m][mkmer]['positive'][name].append(i)
                else:
                    pdict[m][mkmer]['positive'][name] = [i]
    fp.close()
    return targets, pdict


def get_negative_positions(pdict,infile,k,masks):
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

            for m in masks:
                mkmer = ''.join([kmer[ind] for ind in range(k) if ind not in m])


                if mkmer in pdict[m]:
                    if name in pdict[m][mkmer]['negative']:
                        pdict[m][mkmer]['negative'][name].append(i)
                    else:
                        pdict[m][mkmer]['negative'][name] = [i]
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
    masks = [c for i in range(num_degen+1) for c in itertools.combinations(range(plen),i)]
    print("{} masks".format(len(masks)))
    targets, probes_dict = get_all_positions(target_file, plen, masks)

    tot_pos = 0
    for m in masks:
        tot_pos += len(probes_dict[m])
    print("{} positive positions".format(tot_pos))

    # TODO: this assumes that the dictionary for the negative set will be too large,
    # i.e if it is all sequences from a huge metagenome.
    # if that assumption is not correct, then we can just build another dictionary for the negative
    # set instead of doing this.
    get_negative_positions(probes_dict,neg_file,plen,masks)

    tot_neg = 0
    for m in masks:
        for k in probes_dict[m]:
            if len(probes_dict[m][k]['negative'])>0:
                tot_neg += 1
    print("{} negative".format(tot_neg))


    not_allowed = set([(m,k) if len(probes_dict[m][k]['negative'][n]) > num_fp for m in probes_dict for k in probes_dict[m] for n in probes_dict[m][k]['negative']])
    print("{} not allowed:".format(len(not_allowed)))

    # sort the kmers by decreasing number of sequences hit
    kmers_list = [(len(probes_dict[m][k]['positive']),len(probes_dict[m][k]['negative']),m,k) for m in probes_dict for k in probes_dict[m] if (m,k) not in not_allowed]
    kmers_list.sort(key=lambda x: x[1])
    # kmers_list.sort(key=lambda x: x[0],reverse = True) # stable sort - retains order from previous
    # print(kmers_list[:25])
    kmers_list.sort(key=lambda x: x[0]-x[1], reverse=True)
    print(kmers_list[:25])


    # TODO: THIS SHOULD BE A FUNCTION
    covered_targets = set()
    ordered_probes =  []

    hits_dict = {}
    neg_hits_dict = {}

    for kmer in kmers_list:
        pos_positions = probes_dict[kmer[2]][kmer[3]]['positive']
        neg_positions = probes_dict[kmer[2]][kmer[3]]['negative']
        for negative in neg_positions:
            if negative not in neg_hits_dict:
                neg_hits_dict[negative] = neg_positions[negative]
            else:
                neg_hits_dict[negative] = list(merge(neg_hits_dict[negative],neg_positions[negative]))
                if len(neg_hits_dict[negative])>num_fp:
                    # Don't add this probe
                    not_allowed.add(([kmer[2],kmer[3]))
                    break
        if ([kmer[2],kmer[3]) in not allowed: continue # was just added

        # first check if the kmer conflicts with already selected ones
        for target in pos_positions:
            if target not in hits_dict:
                hits_dict[target] = pos_positions[target] # this is a sorted list
            else:
                for p in pos_positions[target]:
                    to_insert = bisect.bisect_left(hits_dict[target],p)
                    if (to_insert>0 and p < hits_dict[target][to_insert-1]+k+MIN_DIST) \
                        or (to_insert<len(hits_dict[target]) and hits_dict[target][to_insert] < p+k+MIN_DIST):
                        # skip inserting this one - count how many times it is inserted, and if it's more than one - add it to ordered_probes
                        not_allowed.add(([kmer[2],kmer[3]))
                    else:
                         hits_dict[target].insert(to_insert,p)
                         if len(hits_dict[target]) > num_probes: covered_targets.add(target)
                         num_hit+=1

        num_hit = 0

        if num_hit > 0: ordered_probes.append(kmer)
        if len(covered_targets) == len(targets): break

    # got here without covering everything - start including false positives
    if len(covered_targets) < len(targets):
        print("{} out of {} targets covered without any FPs".format(len(covered_targets),len(targets)))
        # do the same thing but allow FPs
        kmers_list = [(len(probes_dict[k[2]][k[3]]['positive']),len(probes_dict[k[2]][k[3]]['negative']),k[2],k[3]) for k in not_allowed]
        kmers_list.sort(key=lambda x: len([t for t in hits_dict if t in probes_dict[x[2]][x[3]]['positive']]) # the number of already selected probes it appears on
        # print(kmers_list[:25])
        kmers_list.sort(key=lambda x: x[0]-x[1], reverse=True)
        for kmer in kmers_list:
            pos_positions = probes_dict[kmer[2]][kmer[3]]['positive']
            num_hit = 0
            for target in pos_positions:

                if target not in hits_dict:
                    hits_dict[target] = pos_positions[target] # this is a sorted list
                else:
                    for p in pos_positions[target]:
                        to_insert = bisect.bisect_left(hits_dict[target],p)
                        if (to_insert>0 and p < hits_dict[target][to_insert-1]+k+MIN_DIST) \
                            or (to_insert<len(hits_dict[target]) and hits_dict[target][to_insert] < p+k+MIN_DIST):
                            # skip inserting this one - count how many times it is inserted, and if it's more than one - add it to ordered_probes
                            pass
                        else:
                             hits_dict[target].insert(to_insert,p)
                             if len(hits_dict[target]) > num_probes: covered_targets.add(target)
                             num_hit+=1
            if num_hit > 0: ordered_probes.append(kmer)
            if len(covered_targets) == len(targets): break


    # TODO: DEAL WITH DEGENS







if __name__=='__main__':
    args = parse_user_input()
    main(args)
