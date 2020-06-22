# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 12:49:05 2020

@author: rjovelin
"""

import gzip
import os
from itertools import zip_longest
import regex
import argparse
import json
import time
from barcodex import extract_barcodes