#!/bin/bash


nucmer -t 8 -l 100 --prefix out ankole_WTDBG2.fasta ankole_CANU.fasta
delta-filter -r -q -l 10000 out.delta > out.rq.delta

