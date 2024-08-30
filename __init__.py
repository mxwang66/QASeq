#!/usr/bin/env python3
# -*- coding: utf-8 -*-


__author__ = 'Michael X. Wang'


"""
QASeq analysis code for short/long amplicon panels. 
"""

import os
import sys
CURRENT_DIR = os.path.dirname(__file__)
sys.path.append(CURRENT_DIR)

from QASeq_analysis import main as analysis