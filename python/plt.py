#! /usr/bin/env python3
# pylint: disable=C0103,E0602
"""
...
"""

from array import array

# import this here, otherwise it overwrites the argparse stuff
from ROOT import TColor, TEfficiency, TH1,gStyle, TFile,kOrange, TLine, \
    TCanvas, TLegend, TLatex, TPad
from os import path, makedirs
import cfg
import opt

def get_palette():
    """
    Function to define color palette for 2D efficiency plots
    """

    palette = array('i')

    for i_bin in range(0, 100):
        rgb = []

        if i_bin < 70:
            rgb = [0.70+0.007*i_bin, 0.00+0.0069*i_bin, 0.00]
        elif i_bin < 90:
            rgb = [0.70+0.007*i_bin, 0.00+0.0069 *
                   i_bin+0.10+0.01*(i_bin-70), 0.00]
        else:
            rgb = [0.98-0.098*i_bin-90, 0.80, 0.00]

        palette.append(TColor.GetColor(rgb[0], rgb[1], rgb[2]))

    return palette


class Plotter:
    """
    Handle plotting of one canvas
    """
    canvas = TCanvas('canvas', 'canvas', 800, 800)


    def __init__(self, config):
        """
        Constructor
        """
        self.config = config
        self.options = opt.Options()
        self.options.assign(config.plot.option)

        histos = []
        histo_classes = []

        TH1.AddDirectory(False)  # CB find a better place ...
        gStyle.SetOptTitle(0)
        gStyle.SetPaintTextFormat("1.2f")
        gStyle.SetHistMinimumZero()

        for input_cfg in self.config.inputs:
            file = TFile.Open(input_cfg.file_name)

            if not file:
                raise RuntimeError(f'File not found: {input_cfg.file_name}')

            folder = file.GetDirectory(input_cfg.folder)
            name = None

            for keys in folder.GetListOfKeys():
                if keys.GetName() == input_cfg.plot_name:
                    name = keys.GetName()

                    if name.find("=") > 0 or folder.Get(name).ClassName() == "TDirectoryFile":
                        continue

                    histo = folder.Get(name).Clone()  # CB fix Clone string

                    h_dim = histo.GetDimension()
                    h_class = histo.ClassName()

                    if h_class == "TEfficiency" and h_dim == 1:
                        histo.Paint("")
                        histos.append(histo.GetPaintedGraph())
                    elif h_class == "TEfficiency" and h_dim == 2:
                        histo.Paint("")
                        histos.append(histo.GetPaintedHistogram())
                    else:
                        histos.append(histo)

                    if h_class == "TEfficiency":
                        h_class += f"{histo.GetDimension()}D"
                    histo_classes.append(h_class)

            file.Close()

        self.histo_classes = histo_classes
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

        if axis_z:
            pad.SetRightMargin(0.15)
            self.canvas.SetRightMargin(0.15)

        for i, (histo, h_class) in enumerate(zip(self.histos,self.histo_classes)):

            histo.SetLineWidth(2)
            histo.SetLineColor(config.colors[i])
            histo.SetMarkerStyle(config.markers[i])
            histo.SetMarkerColor(config.colors[i])

            histo.SetTitle(f';{axis_x.label};{axis_y.label};{axis_z.label}')

            if options.chamberSummary:
                if i == 0:
                    dt_color = TColor.GetColorTransparent(kOrange-2, 0.5)
                    histo.SetFillColor(dt_color)
                else:
                    histo.SetLineStyle(7)

            VALID_PLOT_CLASSES = ["TH1I", "TH1D", "TH1F", "TH2F", "TProfile", "TEfficiency1D", "TEfficiency2D"]
            if h_class not in VALID_PLOT_CLASSES:
                raise ValueError(f"h_class: {h_class} not in {VALID_PLOT_CLASSES}")

            order_option = "" if i == 0 else "same"

            if h_class == "TEfficiency1D":
                if i == 0:
                    order_option += "A"
                histo.SetMarkerSize(1.5)
                histo.Draw(f'{options.option}{order_option}P')
                if options.cleanEmptyBins:
                    for i_point in reversed(range(0, histo.GetN())):
                        pointY = histo.GetY()[i_point]
                        if abs(pointY) < 0.005:
                            histo.RemovePoint(i_point)
                    histo.Draw(f'{options.option}{order_option}P')
            else:
                histo.Draw(f'{options.option}{order_option}')

            self.canvas.Update()

            if axis_x:
                histo.GetXaxis().SetRangeUser(axis_x.min, axis_x.max)

            if h_class == "TH1F" or h_class == "TH1D" or h_class == "TH1I":
                if options.scale:
                    histo.Scale(1.0/histo.Integral())
                y_range = [0.0, histo.GetMaximum() * 1.5]
                if axis_y and options.logY:
                    y_range = [axis_y.min, axis_y.max]
                histo.GetYaxis().SetRangeUser(y_range[0], y_range[1])
            elif h_class in ["TProfile", "TEfficiency1D"]:
                histo.GetYaxis().SetRangeUser(axis_y.min, axis_y.max)
            elif h_class == "TEfficiency2D":
                if axis_z:
                    histo.SetMinimum(axis_z.min)
                    histo.SetMaximum(axis_z.max)

                palette = get_palette()
                gStyle.SetPalette(len(palette), palette)
                histo.SetContour(len(palette))

                pad.SetFrameFillColor(851)
                pad.Update()
                histo.GetZaxis().SetLabelSize(0.03)
                histo.Draw(options.option)


                if options.gridBySector:
                    self.canvas.cd(1)
                    line = TLine()
                    for x in range(1, 12):
                        line.DrawLine(x+0.5, -2.5, x+0.5, 2.5)
                    for y in range(-2, 2):
                        line.DrawLine(0.5, y+0.5, 12.5, y+0.5)

                if options.gridByWheel:
                    self.canvas.cd(1)
                    line = TLine()
                    for y in range(1, 4):
                        line.DrawLine((y*5)+1.0, 1.0, (y*5)+1.0, 13.0)
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
                    for y in range(1, 13):
                        histo.GetYaxis().SetBinLabel(y, str(y))

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

            histo.GetYaxis().SetRangeUser(axis_y.min, axis_y.max)

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


            if axis_x.label == 'position' and histo.GetXaxis().GetNbins() <= 21:
                for i_bin in range(1, histo.GetXaxis().GetNbins()+1):
                    histo.GetXaxis().SetBinLabel(i_bin, "WH"+str(((i_bin-1)%5)-2))          #"MB"+str(int((i_bin-1)/5)+1)+
                    histo.GetXaxis().SetTitle(" ")
                    #line = TLine()
                    #for y in range(1, 4):
                    #    x = (y*5)+.3
                    #    line.DrawLine(x , 0.0, x, 3.0)
                    latex = TLatex()
                    latex.SetNDC()
                    latex.SetTextFont(62)
                    latex.SetTextSize(0.028)
                    latex.SetTextAlign(11)
                    latex.DrawLatex(0.15, 0.15, "MB1")
                    latex.DrawLatex(0.35, 0.15, "MB2")
                    latex.DrawLatex(0.55, 0.15, "MB3")
                    latex.DrawLatex(0.75, 0.15, "MB4")
                    histo.GetXaxis().LabelsOption('v')
                    histo.GetXaxis().SetLabelSize(18)
            
            histo.GetYaxis().SetLabelSize(22)
            histo.GetYaxis().SetTitleFont(63)
            histo.GetYaxis().SetLabelFont(43)
            histo.GetYaxis().SetTitleSize(22)
            histo.GetYaxis().SetLabelSize(20)
            histo.GetYaxis().SetTitleOffset(1.5)

            if axis_y.label == 'wheel' and histo.GetYaxis().GetNbins() == 5:
                for i_bin in range(1, histo.GetYaxis().GetNbins()+1):
                    histo.GetYaxis().SetBinLabel(i_bin, str(i_bin - 3))

            if axis_z:
                histo.GetZaxis().SetTitleFont(63)
                histo.GetZaxis().SetTitleSize(22)
                histo.GetZaxis().SetTitleOffset(1.75)

            self.canvas.Update()

        self.canvas.Update()

        self.canvas.cd()

    def __legend__(self):
        """
        Draw legend method
        """

        legend_range = self.config.plot.legend_range
        options = self.options

        self.canvas.cd()

        if len(self.config.inputs) == 1:
            return

        self.canvas.cd(1)
        leg = TLegend(legend_range[0], legend_range[1],
                      legend_range[2], legend_range[3])

        for histo, c in zip(self.histos, self.config.inputs):
            entry = c.name

            if options.fitGaus:
                mean = histo.GetMean()
                rms = histo.GetRMS()
                fit_f = TF1("fitFunc", "gaus", mean - 1.2 * rms, mean + 1.2 * rms)
                histo.Fit("fitFunc", "RMQN")

                if fit_f:
                    mean = fit_f.GetParameter(1)
                    rms = fit_f.GetParameter(2)
                    entry = f"{entry}  [fit mean {mean:3.1f} - #sigma {mean:3.1f} ns]"

            leg.AddEntry(histo, entry, "lep")

        leg.SetBorderSize(1)
        leg.SetTextFont(43)
        leg.SetTextSize(20)

        self.leg = leg
        self.leg.Draw()
        self.canvas.Update()

    def __caption__(self):
        """
        Draw caption method
        """

        config = self.config.plot

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
        Write canvas to file(s)
        """
        folder = self.config.output.directory

        if not path.exists(folder):
            makedirs(folder)

        for file in self.config.output.paths():
            self.canvas.SaveAs(file)


    def Print(self):
        """
        Main function drawing and saving a canvas
        """
        self.__draw__()
        self.__legend__()
        self.__caption__()
        self.__write__()
