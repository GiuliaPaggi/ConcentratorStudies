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
    base_jsons = json.loads(f.read())

    base_cfgs = {}
    for key, item in base_jsons.items():
        print(f'Processing base: {key}')
        base_cfgs[key] = cfg.BaseConfig()
        base_cfgs[key].assign(item)


# Go through plots defined in config JSON

for key, item in cfg_json.items():
    print(f'Processing: {key}')

    plotter_cfg = cfg.Config()
    plotter_cfg.assign(item)
    if plotter_cfg.base not in base_cfgs.keys():
        raise RuntimeError(f'Base configuration: {plotter_cfg.base} not found')
    plotter_cfg += base_cfgs[plotter_cfg.base]
    plotter = plt.Plotter(plotter_cfg)
    plotter.Print()
