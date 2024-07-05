#!/usr/bin/python
import os
import random

def random_barcode(length=18):
    return ''.join(random.choices('ACTG', k=length))

for _ in range(3):
    os.system(f'touch {random_barcode()}.bam')