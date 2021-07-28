# probes
```
Output possible probe sets

optional arguments:
  -h, --help            show this help message and exit
  -t TARGET, --target TARGET
                        fasta file of the target sequences to find probes for
                        (default: None)
  -n NEGATIVE, --negative NEGATIVE
                        fasta file of the negtive sequence set that probes
                        cannot hit (default: None)
  -o OUTDIR, --outdir OUTDIR
                        output directory (default: None)
  -l PROBE_LEN, --probe_len PROBE_LEN
                        length of probes (default: 20)
  -d MAX_DEGENERATE, --max_degenerate MAX_DEGENERATE
                        maximum number of degenerate bases allowed in a probe
                        (default: 1)
  -c MIN_COVERAGE, --min_coverage MIN_COVERAGE
                        minimum number of probes required to consider a target
                        as covered (default: 15)
  -f MAX_FALSE, --max_false MAX_FALSE
                        maximum number of times a negative sequence may be hit
                        (default: 3)
```
