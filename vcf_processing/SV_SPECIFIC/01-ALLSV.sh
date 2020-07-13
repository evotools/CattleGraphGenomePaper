#!/bin/bash

# Extract vcf with abs(SVLEN) > 500 bp
bcftools merge -l samples.txt -O v | sed 's/hereford.//g' | python SVLEN.py - 500 | bgzip -c > JOINT/JOINED.SV500.vcf.gz

