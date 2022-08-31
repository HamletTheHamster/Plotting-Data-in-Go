package main

import (
  "image/color"
  "github.com/Arafatk/glot"
  "github.com/maorshutman/lm"
  "encoding/csv"
  "bufio"
  "fmt"
  "os"
  "io"
  "strconv"
  "strings"
  "math"
  "gonum.org/v1/plot"
  "gonum.org/v1/plot/plotter"
  "gonum.org/v1/plot/vg"
  "gonum.org/v1/plot/font"
  "gonum.org/v1/plot/vg/draw"
  "time"
)

func main() {

  label, file := readMeta()

  pras, pas, prs, ps := getAllData(file, label)
  prasLabel, pasLabel, prsLabel, psLabel := getAllLabels(label)

  raw := []int{}
  if len(raw) > 0 {
    plotRaw(
      raw,
      pras, pas, prs, ps,
      prasLabel, pasLabel, prsLabel, psLabel,
    )
  }

  s, as := subtractBackground(pras, pas, prs, ps)

  subtracted := []int{0,1,2}
  if len(subtracted) > 0 {
    plotSubtracted(subtracted, s, as, prsLabel, prasLabel)
  }

  subtractedTogether := []int{}
  if len(subtractedTogether) > 0 {
    plotSubtractedTogether(subtractedTogether, s, as, prsLabel, prasLabel)
  }

  subtractedGrouped := []int{}
  if len(subtractedGrouped) > 0 {
    goPlotSubGrpd(subtractedGrouped, s, as, prsLabel, prasLabel)
  }

  // Lorentz fit better
  fitSets := []int{0,1,2}
  if len(fitSets) > 0 {

    // Fit parameter guesses
    amp := 2.5
    wid := 0.2
    cen := 2.26

    // as
    fmt.Println("\nAnti-Stokes\n")
    var asFit [][]float64
    var asFits [][][]float64
    var asWidthLine [][]float64
    var asWidthLines [][][]float64
    var asfwhm []float64

    for _, set := range fitSets {

      f := func(dst, guess []float64) {

        amp, wid, cen := guess[0], guess[1], guess[2]

        for i := 0; i < len(as[set][0]); i++ {
          x := as[set][0][i]
          y := as[set][1][i]
          dst[i] = amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)) - y
        }
      }

      jacobian := lm.NumJac{Func: f}

      // Solve for fit
      toBeSolved := lm.LMProblem{
        Dim:        3,
        Size:       len(as[set][0]),
        Func:       f,
        Jac:        jacobian.Jac,
        InitParams: []float64{amp, wid, cen},
        Tau:        1e-6,
        Eps1:       1e-8,
        Eps2:       1e-8,
      }

      results, _ := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

      amp, wid, cen := results.X[0], math.Abs(results.X[1]), results.X[2]

      asfwhm = append(asfwhm, wid*2000)

      fmt.Printf("set %d | width: %.2f MHz | peak: %.6f nV | center: %.4f GHz\n", set, wid*2000, amp, cen)

      var asyFits []float64

      // Create function according to solved fit parameters
      for i := 0; i < len(as[set][0]); i++ {
        // (amp*wid^2/((x-cen)^2+wid^2))
        x := as[set][0][i]
        asyFits = append(asyFits, amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)))
      }

      // Width lines
      asWidthLine = [][]float64{{cen - wid, cen + wid},{amp/2, amp/2}}
      asWidthLines = append(asWidthLines, asWidthLine)

      asFit = [][]float64{as[set][0], asyFits}
      asFits = append(asFits, asFit)
    }

    // goPlot as fits
    goPlotasFits(fitSets, as, prasLabel, asFits, asWidthLines, asfwhm)

    // goPlot power vs width
    goPlotasPowerVsWid(fitSets, prasLabel, asfwhm)

    /*
    // s
    fitStokes := []int{}
    if len(fitStokes) > 0 {

      fmt.Println("\nStokes\n")
      var sFit [][]float64
      var sFits [][][]float64
      var sWidthLine [][]float64
      var sWidthLines [][][]float64
      var sfwhm []float64

      for _, set := range fitSets {

        f := func(dst, guess []float64) {

          amp, wid, cen := guess[0], guess[1], guess[2]

          for i := 0; i < len(s[set][0]); i++ {
            x := s[set][0][i]
            y := s[set][1][i]
            dst[i] = amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)) - y
          }
        }

        jacobian := lm.NumJac{Func: f}

        // Solve for fit
        toBeSolved := lm.LMProblem{
      	  Dim:        3,
       	  Size:       len(s[set][0]),
       	  Func:       f,
       	  Jac:        jacobian.Jac,
       	  InitParams: []float64{amp, wid, cen},
       	  Tau:        1e-6,
       	  Eps1:       1e-8,
       	  Eps2:       1e-8,
        }

        results, _ := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

        amp, wid, cen := results.X[0], math.Abs(results.X[1]), results.X[2]

        sfwhm = append(sfwhm, wid*2000)

        fmt.Printf("set %d | width: %.2f MHz | peak: %.6f nV | center: %.4f GHz\n", set, wid*2000, amp, cen)

        var syFits []float64

        // Create function according to solved fit parameters
        for i := 0; i < len(s[set][0]); i++ {
          // (amp*wid^2/((x-cen)^2+wid^2))
          x := s[set][0][i]
          syFits = append(syFits, results.X[0] * math.Pow(results.X[1], 2) / (math.Pow(x - results.X[2], 2) + math.Pow(results.X[1], 2)))
        }

        // Width lines
        sWidthLine = [][]float64{{cen - wid, cen + wid},{amp/2, amp/2}}
        sWidthLines = append(sWidthLines, sWidthLine)

        sFit = [][]float64{s[set][0], syFits}
        sFits = append(sFits, sFit)
      }
      fmt.Printf("\n")
    }*/
  }
}

func readMeta() ([]string, []string) {

  // Read
  metaFile, err := os.Open("Data/meta.csv")
  if err != nil {
    fmt.Println(err)
  }

  reader := csv.NewReader(metaFile)
  meta, err := reader.ReadAll()
  if err != nil {
    fmt.Println(err)
  }

  var label, data []string

  for _, value := range meta {
    label = append(label, value[1])
    data = append(data, value[3])
  }

  return label, data
}

func getAllData(fileNames, labels []string) ([][][]float64, [][][]float64, [][][]float64, [][][]float64) {

  var pras, pas, prs, ps [][][]float64

  // Assign data by checking csv name
  for _, fileName := range fileNames {
    if strings.Contains(fileName, "ras.") {
      pras = append(pras, getData(&fileName))
    } else if strings.Contains(fileName, "bas.") {
      pas = append(pas, getData(&fileName))
    } else if strings.Contains(fileName, "rs.") {
      prs = append(prs, getData(&fileName))
    } else if strings.Contains(fileName, "bs.") {
      ps = append(ps, getData(&fileName))
    }
  }

  return pras, pas, prs, ps
}

func getData(csvName *string) ([][]float64) {

  // Read
  f, err := os.Open("Data/" + *csvName)
  if err != nil {
    panic(err)
  }
  defer f.Close()
  dataStr, err := readCSV(f)
  if err != nil {
    panic(err)
  }

  // Separate, Strip, & Transpose
  var frequencyStrT, signalStrT []string

  for i := 1; i < 602; i++ {
    frequencyStrT = append(frequencyStrT, strings.ReplaceAll(dataStr[i][0]," ",""))
    signalStrT = append(signalStrT, strings.ReplaceAll(dataStr[i][2]," ",""))
  }

  // Convert to float
  var frequency, signal []float64

  for _, freqElem := range frequencyStrT {
    freqValue, err := strconv.ParseFloat(freqElem, 64)
    if err == nil {
      frequency = append(frequency, freqValue/1e9)
    }
    if err != nil {
      fmt.Println(err)
    }
  }

  for _, sigElem := range signalStrT {
    sigValue, err := strconv.ParseFloat(sigElem, 64)
    if err == nil {
      signal = append(signal, sigValue)
    }
  }

  // Convert to Linear if dBm
  if dataStr[1][3] == " dBm" {
    var nV []float64

    for _, dBm := range signal {
      nV = append(nV, 1000*math.Pow(10, 6)*math.Pow(10, dBm/10.))
    }

    return [][]float64{frequency, nV}
  } else if dataStr[1][3] == "  uV" {
    var nV []float64

    for _, uV := range signal {
      nV = append(nV, 1000*uV)
    }

    /* Conver to picovolts
    var pV []float64
    for _, uV := range signal {
      pV = append(pV, 1000*uV)
    }*/

    return [][]float64{frequency, nV}
  }

  fmt.Println("Warning: check units - not uV or dBm")
  return [][]float64{frequency, signal}
}

func readCSV(rs io.ReadSeeker) ([][]string, error) {
  // Skip first row (line)
  row1, err := bufio.NewReader(rs).ReadSlice('\n')
  if err != nil {
    return nil, err
  }
  _, err = rs.Seek(int64(len(row1)), io.SeekStart)
  if err != nil {
    return nil, err
  }

  // Read remaining rows
  r := csv.NewReader(rs)
  rows, err := r.ReadAll()
  if err != nil {
    return nil, err
  }
  return rows, nil
}

func getAllLabels(label []string) ([]string, []string, []string, []string) {

  var prasLabel, pasLabel, prsLabel, psLabel []string

  // Assign labels by checking label
  for _, thisLabel := range label {
    if strings.Contains(thisLabel, "pras") {
      prasLabel = append(prasLabel, thisLabel)
    } else if strings.Contains(thisLabel, "pas") {
      pasLabel = append(pasLabel, thisLabel)
    } else if strings.Contains(thisLabel, "prs") {
      prsLabel = append(prsLabel, thisLabel)
    } else if strings.Contains(thisLabel, "ps") {
      psLabel = append(psLabel, thisLabel)
    }
  }

  return prasLabel, pasLabel, prsLabel, psLabel
}

func buildData(data [][]float64) (plotter.XYs) {

  xy := make(plotter.XYs, len(data[0]))

  for i := range xy {
    xy[i].X = data[0][i]
    xy[i].Y = data[1][i]
  }

  return xy
}

func plotRaw(
  sets []int,
  pras, pas, prs, ps [][][]float64,
  prasLabel, pasLabel, prsLabel, psLabel []string,
  ) {

  for i := 0; i < len(sets); i++ {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Raw")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    plot.AddPointGroup(prasLabel[sets[i]], "points", pras[sets[i]])
    plot.AddPointGroup(pasLabel[sets[i]], "points", pas[sets[i]])
    plot.AddPointGroup(prsLabel[sets[i]], "points", prs[sets[i]])
    plot.AddPointGroup(psLabel[sets[i]], "points", ps[sets[i]])
  }
}

func subtractBackground(pras, pas, prs, ps [][][]float64) ([][][]float64, [][][]float64) {

  var s, as [][][]float64

  for i := 0; i < len(pras); i++ {
    s = append(s, subtract(ps[i], prs[i]))
    as = append(as, subtract(pas[i], pras[i]))
  }

  return s, as
}

func subtract(b, s [][]float64) ([][]float64) {

  var sum float64
  n := 100

  for i := 0; i < n; i++ {
    sum += b[1][i] - s[1][i]
  }

  for i := 0; i < len(b[0]); i++ {
    s[1][i] = s[1][i] - b[1][i] + sum/float64(n)
  }

  return s
}

func plotSubtracted(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  for _, set := range sets {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Background Subtracted")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (nV)")

    plot.AddPointGroup(strings.Trim(sLabel[set], " rs") + " s", "points", s[set])
    plot.AddPointGroup(strings.Trim(asLabel[set], " ras") + " as", "points", as[set])
  }
}

func plotSubtractedTogether(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (nV)")

  for _, set := range sets {
    plot.AddPointGroup(strings.Trim(sLabel[set], " prs") + " s", "points", s[set])
    plot.AddPointGroup(strings.Trim(asLabel[set], " pras") + " as", "points", as[set])
  }
}

func goPlotSubGrpd(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  // as
  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Anti-Stokes Probe Spectra"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 2
  p.X.Max = 2.5
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 36
  p.X.Tick.Label.Font.Variant = "Sans"

  p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 2.05, Label: ""},
    {Value: 2.1, Label: "2.1"},
    {Value: 2.15, Label: ""},
    {Value: 2.2, Label: "2.2"},
    {Value: 2.25, Label: ""},
    {Value: 2.3, Label: "2.3"},
    {Value: 2.35, Label: ""},
    {Value: 2.4, Label: "2.4"},
    {Value: 2.45, Label: ""},
    {Value: 2.5, Label: "2.5"},
  })
  p.X.Padding = vg.Points(-12.5)

  p.Y.Label.Text = "Spectral Density (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = -1
  p.Y.Max = 2.75
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: -.75, Label: ""},
    {Value: -.5, Label: "-.5"},
    {Value: -.25, Label: ""},
    {Value: 0, Label: "0"},
    {Value: .25, Label: ""},
    {Value: .5, Label: ".5"},
    {Value: .75, Label: ""},
    {Value: 1, Label: "1"},
    {Value: 1.25, Label: ""},
    {Value: 1.5, Label: "1.5"},
    {Value: 1.75, Label: ""},
    {Value: 2, Label: "2"},
    {Value: 2.25, Label: ""},
    {Value: 2.5, Label: "2.5"},
    {Value: 2.75, Label: ""},
    {Value: 3, Label: "3"},
  })
  p.Y.Padding = vg.Points(-4.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)
  p.Legend.Add("Pump")

  setColors := make([]color.RGBA, len(sets))
  setColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
  setColors[1] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
  setColors[2] = color.RGBA{R: 122, G: 156, B: 255, A: 255}

  setFitColors := make([]color.RGBA, len(sets))
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[1] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[2] = color.RGBA{R: 99, G: 124, B: 198, A: 255}


  for _, set := range sets {

    asPts := buildData(as[set])

    // Make a scatter plotter and set its style.
    plotSet, err := plotter.NewScatter(asPts)
    if err != nil {
      panic(err)
    }

    plotSet.GlyphStyle.Color = setColors[set]
    plotSet.GlyphStyle.Radius = vg.Points(3)
    plotSet.Shape = draw.CircleGlyph{}

    p.Add(plotSet)

    // Legend
    l, err := plotter.NewScatter(asPts)
    if err != nil {
      panic(err)
    }

    l.GlyphStyle.Color = setFitColors[set]
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}
    p.Legend.Add(strings.Trim(asLabel[set], " ras"), l)
  }

  date := time.Now()

  // Make current date folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02"), 0755); err != nil {
      panic(err)
    }
  }

  // Make current time folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05"), 0755); err != nil {
      panic(err)
    }
  }

  savePlotAs := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/Anti-Stokes Background Subtracted"
  // Save the plot to a PNG file.
  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".png"); err != nil {
  	panic(err)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".svg"); err != nil {
    panic(err)
  }
}

func goPlotasFits(sets []int, as [][][]float64, labels []string,
  fits [][][]float64, widthLines [][][]float64, widths []float64) {

    p := plot.New()
    p.BackgroundColor = color.RGBA{A:0}
    p.Title.Text = "Anti-Stokes Probe Spectra"
    p.Title.TextStyle.Font.Typeface = "liberation"
    p.Title.TextStyle.Font.Variant = "Sans"
    p.Title.TextStyle.Font.Size = 50
    p.Title.Padding = font.Length(50)

    p.X.Label.Text = "Frequency (GHz)"
    p.X.Label.TextStyle.Font.Variant = "Sans"
    p.X.Label.TextStyle.Font.Size = 36
    p.X.Label.Padding = font.Length(20)
    p.X.LineStyle.Width = vg.Points(1.5)
    p.X.Min = 2
    p.X.Max = 2.5
    p.X.Tick.LineStyle.Width = vg.Points(1.5)
    p.X.Tick.Label.Font.Size = 36
    p.X.Tick.Label.Font.Variant = "Sans"

    p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: 2.05, Label: ""},
      {Value: 2.1, Label: "2.1"},
      {Value: 2.15, Label: ""},
      {Value: 2.2, Label: "2.2"},
      {Value: 2.25, Label: ""},
      {Value: 2.3, Label: "2.3"},
      {Value: 2.35, Label: ""},
      {Value: 2.4, Label: "2.4"},
      {Value: 2.45, Label: ""},
      {Value: 2.5, Label: "2.5"},
    })
    p.X.Padding = vg.Points(-12.5)

    p.Y.Label.Text = "Spectral Density (nV)"
    p.Y.Label.TextStyle.Font.Variant = "Sans"
    p.Y.Label.TextStyle.Font.Size = 36
    p.Y.Label.Padding = font.Length(20)
    p.Y.LineStyle.Width = vg.Points(1.5)
    p.Y.Min = -1
    p.Y.Max = 2.75
    p.Y.Tick.LineStyle.Width = vg.Points(1.5)
    p.Y.Tick.Label.Font.Size = 36
    p.Y.Tick.Label.Font.Variant = "Sans"
    p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: -.75, Label: ""},
      {Value: -.5, Label: "-.5"},
      {Value: -.25, Label: ""},
      {Value: 0, Label: "0"},
      {Value: .25, Label: ""},
      {Value: .5, Label: ".5"},
      {Value: .75, Label: ""},
      {Value: 1, Label: "1"},
      {Value: 1.25, Label: ""},
      {Value: 1.5, Label: "1.5"},
      {Value: 1.75, Label: ""},
      {Value: 2, Label: "2"},
      {Value: 2.25, Label: ""},
      {Value: 2.5, Label: "2.5"},
      {Value: 2.75, Label: ""},
      {Value: 3, Label: "3"},
    })
    p.Y.Padding = vg.Points(-4.75)

    p.Legend.TextStyle.Font.Size = 36
    p.Legend.TextStyle.Font.Variant = "Sans"
    p.Legend.Top = true
    p.Legend.XOffs = vg.Points(-50)
    p.Legend.YOffs = vg.Points(-50)
    p.Legend.Padding = vg.Points(10)
    p.Legend.ThumbnailWidth = vg.Points(50)
    p.Legend.Add("Pump")

    setPtColors := make([]color.RGBA, len(sets))
    setPtColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
    setPtColors[1] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
    setPtColors[2] = color.RGBA{R: 122, G: 156, B: 255, A: 255}

    setFitColors := make([]color.RGBA, len(sets))
    setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    setFitColors[1] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    setFitColors[2] = color.RGBA{R: 99, G: 124, B: 198, A: 255}


    for _, set := range sets {

      pts := buildData(as[set])
      fit := buildData(fits[set])
      wid := buildData(widthLines[set])

      // Plot points
      plotPts, err := plotter.NewScatter(pts)
      if err != nil {
        panic(err)
      }

      plotPts.GlyphStyle.Color = setPtColors[set]
      plotPts.GlyphStyle.Radius = vg.Points(3)
      plotPts.Shape = draw.CircleGlyph{}

      // Plot fit
      plotFit, err := plotter.NewLine(fit)
      if err != nil {
        panic(err)
      }

      plotFit.LineStyle.Color = setFitColors[set]
      plotFit.LineStyle.Width = vg.Points(3)

      // Width lines
      plotWid, err := plotter.NewLine(wid)
      if err != nil {
        panic(err)
      }

      plotWid.LineStyle.Color = setFitColors[set]
      plotWid.LineStyle.Width = vg.Points(4)
      plotWid.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

      // Add set plots to p
      p.Add(plotPts, plotFit, plotWid)

      // Legend
      l, err := plotter.NewScatter(pts)
      if err != nil {
        panic(err)
      }

      l.GlyphStyle.Color = setFitColors[set]
      l.GlyphStyle.Radius = vg.Points(6)
      l.Shape = draw.CircleGlyph{}
      p.Legend.Add(strings.Trim(labels[set], " pras"), l)
    }

    // Save plot
    name := "Anti-Stokes w Fits"
    date := time.Now()

    // Make current date folder if it doesn't already exist
    if _, err := os.Stat("plots/" + date.Format("2006-Jan-02")); os.IsNotExist(err) {
      if err := os.Mkdir("plots/" + date.Format("2006-Jan-02"), 0755); err != nil {
        panic(err)
      }
    }

    // Make current time folder if it doesn't already exist
    if _, err := os.Stat("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")); os.IsNotExist(err) {
      if err := os.Mkdir("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05"), 0755); err != nil {
        panic(err)
      }
    }

    path := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/" + name
    // Save the plot to a PNG file.
    if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".png"); err != nil {
    	panic(err)
    }

    if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".svg"); err != nil {
      panic(err)
    }
}

func goPlotasPowerVsWid(sets []int, labels []string, widths []float64) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Pump Power vs Widths of Fits"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Pump Power (mW)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 0
  p.X.Max = 200
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 36
  p.X.Tick.Label.Font.Variant = "Sans"

  p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 0, Label: "0"},
    {Value: 25, Label: ""},
    {Value: 50, Label: "50"},
    {Value: 75, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 125, Label: ""},
    {Value: 150, Label: "150"},
    {Value: 175, Label: ""},
    {Value: 200, Label: "200"},
  })
  p.X.Padding = vg.Points(-8.25)

  p.Y.Label.Text = "Full Width Half Max (MHz)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 90
  p.Y.Max = 130
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 90, Label: "90"},
    {Value: 95, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 105, Label: ""},
    {Value: 110, Label: "110"},
    {Value: 115, Label: ""},
    {Value: 120, Label: "120"},
    {Value: 125, Label: ""},
    {Value: 130, Label: "130"},
  })
  p.Y.Padding = vg.Points(1)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setFitColors := make([]color.RGBA, len(sets))
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[1] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[2] = color.RGBA{R: 99, G: 124, B: 198, A: 255}

  for _, set := range sets {

    pts := make(plotter.XYs, 1)

    if pow, err := strconv.ParseFloat(strings.Trim(labels[set], " mW pras"), 64); err == nil {
      pts[0].X = pow
    } else {
      panic(err)
    }

    pts[0].Y = widths[set]

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      panic(err)
    }

    plotPts.GlyphStyle.Color = setFitColors[set]
    plotPts.GlyphStyle.Radius = vg.Points(6)
    plotPts.Shape = draw.CircleGlyph{}

    // Dashed eye guide lines
    v := make(plotter.XYs, 2)
    h := make(plotter.XYs, 2)

    // Vertical
    v[0].X = pts[0].X
    v[0].Y = 90
    v[1].X = pts[0].X
    v[1].Y = pts[0].Y

    vDash, err := plotter.NewLine(v)
    if err != nil {
      panic(err)
    }

    vDash.LineStyle.Color = setFitColors[set]
    vDash.LineStyle.Width = vg.Points(4)
    vDash.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

    // Horizontal
    h[0].X = -15
    h[0].Y = pts[0].Y
    h[1].X = pts[0].X
    h[1].Y = pts[0].Y

    hDash, err := plotter.NewLine(h)
    if err != nil {
      panic(err)
    }

    hDash.LineStyle.Color = color.RGBA{R: 127, G: 127, B: 127, A: 255}
    hDash.LineStyle.Width = vg.Points(1)
    hDash.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}

    // Add set plots to p
    p.Add(plotPts, vDash, hDash)
    p.Legend.Add(strings.Trim(labels[set], " pras"), plotPts)
  }

  // Save plot
  name := "Pow vs Wid"
  date := time.Now()

  // Make current date folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02"), 0755); err != nil {
      panic(err)
    }
  }

  // Make current time folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05"), 0755); err != nil {
      panic(err)
    }
  }

  path := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/" + name
  // Save the plot to a PNG file.
  if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".png"); err != nil {
    panic(err)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".svg"); err != nil {
    panic(err)
  }
}

func normalizeFit(fit []float64) ([]float64) {

  var shift float64 = (fit[0] + fit[599])/2

  for i := 0; i < 600; i++ {
    fit[i] = fit[i] - shift
  }
  return fit
}
