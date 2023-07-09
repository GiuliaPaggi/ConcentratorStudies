#!/usr/bin/env python3
"""
This program takes an input JSON config and extracts plots from ROOT files.
The output consists of a plot with superimposed plots from multiple files.
"""

import argparse
import json

from ROOT import gROOT

import plt
import cfg


# set ROOT to batch mode, this suppresses printing canvases
gROOT.SetBatch(True)
# suppress stdout pollution of canvas.Print(...)
gROOT.ProcessLine("gErrorIgnoreLevel = 1001;")

# Setup argument parser
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('configJSON',
                    help='Path to the input JSON config file')
parser.add_argument('-b',
                    '--baseJSON',
                    default='config/configPlotterBase.json',
                    help='Path to the input JSON base config file')
args = parser.parse_args()

# Parse JSON file

with open(args.configJSON, 'r', encoding='utf-8') as f:
    cfg_json = json.loads(f.read())

with open(args.baseJSON,'r', encoding='utf-8') as f:
    base_json = json.loads(f.read())
    base_cfg = cfg.BaseConfig()
    base_cfg.assign(base_json['base'])


# Go through plots defined in config JSON

for key, item in cfg_json.items():
    print(f'Processing: {key}')

    plotter_cfg = cfg.Config()
    plotter_cfg.assign(item)
    plotter_cfg += base_cfg

    plotter = plt.Plotter(plotter_cfg)
    plotter.Print()
