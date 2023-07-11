#! /usr/bin/env python3
# pylint: disable=C0103,E0602
"""
...
"""

import sys
from array import array

from ROOT import TColor, TH1, gStyle, TFile, kOrange, TLine, TCanvas, TLegend, TLatex, TPad
from os import path, makedirs
import cfg
import opt

def add_palette(cls):
    """
    Function to define color palette for 2D efficiency plots
    """
    cls.__eff_palette__ = []

    for i_bin in range(0, 100):
        rgb = []

        if i_bin < 70:
            rgb = [0.70+0.007*i_bin, 0.00+0.0069*i_bin, 0.00]
        elif i_bin < 90:
            rgb = [0.70+0.007*i_bin, 0.00+0.0069 *
                   i_bin+0.10+0.01*(i_bin-70), 0.00]
        else:
            rgb = [0.98-0.098*i_bin-90, 0.80, 0.00]

    cls.__eff_palette__.append(TColor.GetColor(rgb[0], rgb[1], rgb[2]))

    return cls


@add_palette
class Plotter:
    """
    Handle plotting of one canvas
    """
    canvas = TCanvas('canvas', 'canvas', 800, 800)
    config = cfg.Config()
    options = opt.Options()


    def __init__(self, config):
        """
        Constructor
        """
        self.config = config
        self.options.assign(config.plot.option)

        histos = []

        TH1.AddDirectory(False)  # CB find a better place ...
        gStyle.SetOptTitle(0)
        gStyle.SetPaintTextFormat("1.3f")
        gStyle.SetHistMinimumZero()

        for input_cfg in self.config.inputs:
            file = TFile.Open(input_cfg.file_name)

            if not file:  # CB improve
                print(f'[ERROR] File not found: {input_cfg.file_name}')
                sys.exit()
            folder = file.GetDirectory(input_cfg.folder)
            name = None

            for keys in folder.GetListOfKeys():
                if keys.GetName() == input_cfg.plot_name:
                    name = keys.GetName()

                    if name.find("=") > 0 or folder.Get(name).ClassName() == "TDirectoryFile":
                        continue

                    histo = folder.Get(name).Clone()  # CB fix Clone string
                    histos.append(histo)

            file.Close()
        self.histos = histos


    def __draw__(self):
        """
        Draw method
        """
        config = self.config.plot
        options = self.options

        # Setup canvas with all elements
        self.canvas.cd()

        pad = TPad('pad', 'pad', 0.01, 0.00, 1.00, 1.00)

        pad.SetGrid()
        pad.Draw()
        pad.cd()


        axis_x = config.x
        axis_y = config.y
        axis_z = config.z
        
        if options.chamberSummary:
            gStyle.SetOptStat("emu")
            gStyle.SetStatX(0.9)
            gStyle.SetStatY(0.9)
        
        gStyle.SetOptStat(0)

        for i, histo in enumerate(self.histos):

            histo.SetLineWidth(2)
            histo.SetLineColor(config.colors[i])
            histo.SetMarkerStyle(config.markers[i])
            histo.SetMarkerColor(config.colors[i])

            h_dim = histo.GetDimension()

            histo.SetTitle(f';{axis_x.label};{axis_y.label}')

            if options.chamberSummary:
                dt_color = TColor.GetColorTransparent(kOrange-2, 0.5)
                histo.SetFillColor(dt_color)

            if i == 0:
                if h_dim == 2:  # CB improve
                    histo.Draw(options.option)
                elif options.cleanEmptyBins:
                    histo.Paint(f'{options.option}AP')
                    histo = histo.GetPaintedGraph()
                    for i_point in reversed(range(0, histo.GetN())):
                        pointY = histo.GetY()[i_point]
                        if pointY < 0.05:
                            histo.RemovePoint(i_point)
                        histo.Draw(f'{options.option}AP')
                else:
                    histo.Draw(options.option)
            else:
                if options.cleanEmptyBins:
                    histo.Paint(f'same{options.option}')
                    histo = histo.GetPaintedGraph()
                    for i_point in reversed(range(0, histo.GetN())):
                        point_y = histo.GetY()[i_point]
                        if point_y < 0.005:
                            histo.RemovePoint(i_point)
                        histo.Draw(f'same{options.option}P')
                else:
                    histo.Draw(f'same{options.option}')

            self.canvas.Update()

            h_class = histo.ClassName()

            if h_class == "TEfficiency" and h_dim == 1:
                histo = histo.GetPaintedGraph()
            elif h_class == "TEfficiency" and h_dim == 2:
                histo = histo.GetPaintedHistogram()

            histo.GetXaxis().SetRangeUser(axis_x.min, axis_x.max)

            if options.scale and h_class == "TH1F":
                histo.Scale(1.0/histo.Integral())

            if h_class == "TH1F":
                y_range = [axis_y.min, axis_y.max] if options.logY else [0.0, histo.GetMaximum() * 1.5]
                histo.GetYaxis().SetRangeUser(y_range[0], y_range[1])

            if h_class in ["TProfile", "TEfficiency"] and h_dim == 1:
                histo.GetYaxis().SetRangeUser(axis_y.min, axis_y.max)
            elif h_class == "TEfficiency" and h_dim == 2:
                nBins = len(self.__eff_palette__)
                gStyle.SetPalette(nBins, array('i', self.__eff_palette__))

                histo.SetMinimum(axis_z.min)
                histo.SetMaximum(axis_z.max)
                histo.SetContour(nBins)
                histo.Draw(options.option)

                if options.gridBySector:
                    line = TLine()
                    for x in range(1, 12):
                        line.DrawLine(x+0.5, -2.5, x+0.5, 2.5)
                    for y in range(-2, 2):
                        line.DrawLine(0.5, y+0.5, 12.5, y+0.5)

                if options.gridByWheel:
                    line = TLine()
                    for x in range(1, 4):
                        line.DrawLine((x*5)+1.0, 1.0, (x*5)+1.0, 13.0)
                    for x in range(0, 20):
                        histo.GetXaxis().SetBinLabel(x+1, 'YB'+str(x % 5 - 2))
                    histo.GetXaxis().LabelsOption('v')
                    latex = TLatex()
                    latex.SetNDC()
                    latex.SetTextFont(62)
                    latex.SetTextSize(0.028)
                    latex.SetTextAlign(11)
                    latex.DrawLatex(0.19, 0.01, "MB1")
                    latex.DrawLatex(0.39, 0.01, "MB2")
                    latex.DrawLatex(0.59, 0.01, "MB3")
                    latex.DrawLatex(0.79, 0.01, "MB4")
                    histo.GetYaxis().SetNdivisions(histo.GetNbinsY(), True)
                    for x in range(1, 13):
                        histo.GetYaxis().SetBinLabel(x, str(x))

                if options.grid or options.gridBySector or options.gridByWheel:
                    histoClone = histo.Clone()
                    histoClone.GetXaxis().SetNdivisions(histo.GetNbinsX(), True)
                    histoClone.GetYaxis().SetNdivisions(histo.GetNbinsY(), True)
                    histoClone.GetXaxis().SetTickLength(0.01)
                    histoClone.GetYaxis().SetTickLength(0.01)
                    histoClone.Draw('same axig')
            elif h_class == 'TH2F':
                gStyle.SetPalette(1)
                histo.Draw(options.option)
            else:
                histo.GetXaxis().SetTitle(axis_x.label)
                histo.GetYaxis().SetTitle(axis_y.label)

            self.canvas.Update()

            histo.GetXaxis().SetLabelSize(22)
            histo.GetXaxis().SetTitleFont(63)
            histo.GetXaxis().SetLabelFont(43)
            histo.GetXaxis().SetTitleSize(22)
            histo.GetXaxis().SetLabelSize(20)
            histo.GetXaxis().SetTitleOffset(1.2)

            if axis_x.label == 'sector' and histo.GetXaxis().GetNbins() <= 14:
                for i_bin in range(1, histo.GetXaxis().GetNbins()+1):
                    histo.GetXaxis().SetBinLabel(i_bin, str(i_bin))

            histo.GetYaxis().SetLabelSize(22)
            histo.GetYaxis().SetTitleFont(63)
            histo.GetYaxis().SetLabelFont(43)
            histo.GetYaxis().SetTitleSize(22)
            histo.GetYaxis().SetLabelSize(20)
            histo.GetYaxis().SetTitleOffset(1.5)

            if axis_y.label == 'wheel' and histo.GetYaxis().GetNbins() == 5:
                for i_bin in range(1, histo.GetYaxis().GetNbins()+1):
                    histo.GetYaxis().SetBinLabel(i_bin, str(i_bin - 3))

            self.canvas.Update()

        self.canvas.Update()

        self.canvas.cd()

    def __draw_legend__(self):
        """
        Draw legend method
        """

        config = self.config.plot
        options = self.options

        if len(self.histos) > 1:

            leg = TLegend(config.legend_range[0], config.legend_range[1],
                          config.legend_range[2], config.legend_range[3])

            for h, c in zip(self.histos, self.config.inputs):
                entry = c.legend_entry

                if options.fitGaus:
                    mean = h.GetMean()
                    rms = h.GetRMS()
                    fit_f = TF1('fitFunc', 'gaus', mean -
                                1.2 * rms, mean + 1.2 * rms)
                    h.Fit("fitFunc", "RMQN")

                    if fit_f:
                        mean = fit_f.GetParameter(1)
                        rms = fit_f.GetParameter(2)
                        entry = f'{entry}  [fit mean {mean:3.1f} - #sigma {mean:3.1f} ns]'

                leg.AddEntry(entry, entry, 'LP')

            leg.SetBorderSize(1)
            leg.SetTextFont(43)
            leg.SetTextSize(20)
            leg.Draw()

        self.canvas.cd()
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextFont(61)
        latex.SetTextSize(0.030)
        latex.SetTextAlign(11)
        latex.DrawLatex(0.115, 0.91, config.logo[0])
        latex.SetTextFont(52)
        latex.SetTextSize(0.027)
        latex.SetTextAlign(11)
        latex.DrawLatex(0.179, 0.91, config.logo[1])
        latex.SetTextFont(42)
        latex.SetTextSize(0.030)
        latex.SetTextAlign(31)
        latex.DrawLatex(0.90, 0.91, config.caption)
        latex.SetTextAlign(11)
        latex.SetTextColor(1)
        latex.SetTextFont(61)
        latex.SetTextSize(0.032)
        latex.DrawLatex(config.legend_range[0],
                        config.legend_range[3] + 0.02,
                        config.legend_title)
        self.canvas.Update()

    def __write__(self):
        """
        Write to file
        """
        folder = self.config.output.directory

        if not path.exists(folder):
            makedirs(folder)

        for file in self.config.output.paths():
            self.canvas.SaveAs(file)

    def Print(self):
        """
        Man function drawing and saving a canvas
        """
        self.__draw__()
        self.__draw_legend__()
        self.__write__()
