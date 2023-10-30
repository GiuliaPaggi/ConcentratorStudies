#! /usr/bin/env python3
# pylint: disable=C0103
"""
Provides a class handling extra drawing option used by plotter Objects
"Basic" unit-tests:
    are provided as __main__ script
TO-DO:
    - improve assign
"""

from dataclasses import dataclass


@dataclass
class Options :
    """
    Handle display options
    """
    gridByWheel: bool = False
    gridBySector: bool = False
    grid: bool = False
    chamberSummary: bool = False
    cleanEmptyBins: bool = False
    scale: bool = False
    fitGaus: bool = False
    logY: bool = False

    option: str = ''

    def assign(self, opt):
        """
        Assign values
        """
        for name in vars(self).keys():
            if name in opt:
                setattr(self, name, True)
                opt = opt.replace(name, "")

        self.option = opt


if __name__ == '__main__':
    option = "gridByWheelcleanEmptyBinstext"

    options = Options()
    options.assign(option)
    print(options)
