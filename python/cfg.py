#! /usr/bin/env python3
# pylint: disable=C0103
"""
Provides a set of classes handling the configuration of plotter Objects
__add__ operators:
    are overloaded to allow updating 'partial' configurations
    by summing a 'partial' configuration with 'base' configuration.
"Basic" unit-tests:
    are provided as __main__ script
TO-DO:
    - improve unit-tests
    - add an is_complete member function
    - remove redundant cfg parameters
    - improve code
"""

import os
from dataclasses import dataclass, field


@dataclass
class Axis:
    """
    Handle axis information
    """
    min: float = 0.0
    max: float = 0.0
    label: str = ''

    def __nonzero__(self):
        return self.min < self.max

    def assign(self, vec):
        """
        Assign values
        """
        self.min = vec[0]
        self.max = vec[1]
        self.label = vec[2]

    def __str__(self):
        """
        Print
        """
        result = "Axis: {"
        for name, item in vars(self).items():
            result = f"{result} {name}: {item} "
        result = f"{result}}}"
        return result

    def __bool__(self):
        """
        Cast to bool defining a non valid Axis
        """
        return not (self.min == 0.0 and self.max == 0.0 and self.label == "")


# pylint: disable=R0902
@dataclass
class Plot:
    """
    Handle plot information
    """
    x: Axis = field(default_factory=Axis)
    y: Axis = field(default_factory=Axis)
    z: Axis = field(default_factory=Axis)
    option: str = ''
    caption: str = ''
    logo: list[str] = field(default_factory=list)
    colors: list[int] = field(default_factory=list)
    markers: list[int] = field(default_factory=list)
    legend_title: str = ''
    legend_range: list[float] = field(default_factory=list)

    def assign(self, cfg):
        """
        Assign values
        """
        for name, item in vars(self).items():
            if name in cfg.keys():
                if name in ['x', 'y', 'z']:
                    item.assign(cfg[name])
                else:
                    setattr(self, name, cfg[name])

    def __str__(self):
        """
        Print
        """
        result = "Plot: {\n"
        for name, item in vars(self).items():
            result = f"{result}  {name}: {item}\n"
        result = f"{result}}}"
        return result

# pylint: enable=R0902


@dataclass
class Output:
    """
    Handle output information
    """
    directory: str = ''
    file_name: str = ''
    file_types: list[str] = field(default_factory=list)

    def assign(self, cfg):
        """
        Assign values
        """
        for name in vars(self).keys():
            if name in cfg.keys():
                setattr(self, name, cfg[name])

    def paths(self):
        """
        Return list of paths
        """
        return [os.path.join(self.directory, f'{self.file_name}.{t}') for t in self.file_types]

    def __str__(self):
        """
        Print
        """
        result = "Output:\n"
        for name, item in vars(self).items():
            result = f"{result} {name}: {item}\n"
        return result


@dataclass
class Input:
    """
    Handle plot information
    """
    name: str = ''
    file_name: str = 'results/outputFile_DoubleMuon_FlatPt-1To100_noPU_noRPC.root'
    folder: str = ''
    plot_name: str = ''

    def assign(self, n, cfg):
        """
        Assign values
        """
        self.name = n
        for name in vars(self).keys():
            if name in cfg.keys():
                setattr(self, name, cfg[name])

    def __str__(self):
        """
        Print
        """
        result = "Input:\n"
        for item, val in vars(self).items():
            result = f"{result} {item}: {val}\n"
        return result


@dataclass
class BaseConfig:
    """
    Handle base config
    """
    plot: Plot = field(default_factory=Plot)
    output: Output = field(default_factory=Output)

    def assign(self, cfg):
        """
        Assign values
        """
        for name, var in vars(self).items():
            if name in cfg.keys():
                if name in ['plot', 'output']:
                    var.assign(cfg[name])
                else:
                    setattr(self, name, cfg[name])

    def __str__(self):
        """
        Print
        """
        result = "BaseConfig:\n"
        for item, val in vars(self).items():
            result = f"{result} {item}: {val}\n"
        return result


@dataclass
class Config(BaseConfig):
    """
    Handle config
    """
    base: str = ''
    inputs: list[Input] = field(default_factory=list)

    def assign(self, cfg):
        """
        Assign values
        """
        for name, var in vars(self).items():
            if name in cfg.keys():
                if name in ['plot', 'output']:
                    var.assign(cfg[name])
                elif name == 'inputs':
                    for key, i_cfg in cfg[name].items():
                        tmp = Input()
                        tmp.assign(key,i_cfg)
                        getattr(self, name).append(tmp)
                else:
                    setattr(self, name, cfg[name])

    def __add__(self, other):
        """
        Add
        """
        cfg = self

        for category in ['plot', 'output']:
            for name in vars(getattr(cfg, category)).keys():
                if not getattr(getattr(cfg, category), name):
                    attr = getattr(getattr(other, category), name)
                    setattr(getattr(cfg, category), name, attr)

        return cfg

    def __radd__(self, other):
        """
        Add
        """
        if other == 0:
            return self

        return self.__add__(other)

    def __str__(self):
        """
        Print
        """
        result = "Config:\n"
        for item, val in vars(self).items():
            result = f"{result} {item}: {val}\n"
        return result


if __name__ == '__main__':
    a = Axis()
    a.assign([0, 1, 'axis title'])
    print(a)

    plot_cfg = {
        "x": [0.0, 200.0, "mu^{+} mu^{-} inv mass (GeV)"],
        "y": [0.8, 1.05, "Entries"],
        "option": "",
        "legend_title": "leg",
        "caption": "RunD (2022 pp data, 13 TeV)",
        "logo": ["CMS", "Preliminary"],
        "colors": [1, 2, 3, 4, 5, 6, 7, 8, 9],
        "markers": [20, 21, 22, 23, 24, 25, 26],
        "legend_range": [0.49, 0.67, 0.85, 0.80]
    }
    p = Plot()
    p.assign(plot_cfg)
    print(p)

    output_cfg = {
        "directory": "plots/",
        "file_name": "pair_mass",
        "file_types": ["png", "pdf"]
    }

    o = Output()
    o.assign(output_cfg)
    print(o, o.paths())

    input_cfg = {
        "legend_entry": "inv. mass",
        "file_name": "results.root",
        "folder": "folder",
        "plot": "pairMass"
    }

    i = Input()
    i.assign("input",input_cfg)
    print(i)

    base_cfg = {
        'plot': plot_cfg,
        'output': output_cfg
    }

    b = BaseConfig()
    b.assign(base_cfg)
    print(b)

    full_cfg = {
        'plot': plot_cfg,
        'output': output_cfg,
        'inputs': {
            'one': input_cfg,
            'two': input_cfg
        }
    }

    c = Config()
    c.assign(full_cfg)
    print(c)

    full_cfg['plot']['z'] = [0.0, 100.0, "arb. units"]

    c_two = Config()
    c_two.assign(full_cfg)

    base_cfg['plot']['option'] = 'P'

    b_two = BaseConfig()
    b_two.assign(base_cfg)

    c_two = b_two + c_two

    print(c_two)
