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
  "flag"
)

func main() {

  // Flags
  var temp, lcof bool
  flag.BoolVar(&temp, "t", false, "parse sample temps")
  flag.BoolVar(&lcof, "l", false, "liquid-core optical fiber sample")
  flag.Parse()

  if temp {
    fmt.Printf("\n")
    fmt.Println("*Temperature-dependent data*")
  }
  if lcof {
    fmt.Printf("\n")
    fmt.Println("*liquid-core optical fiber sample*")
  }

  date, run, label, file, asNotes, sNotes := readMeta(temp)

  fmt.Printf("\n" + date + " run " + run + "\n")

  ras, bas, rs, bs := getAllData(file, label)
  rasLabel, basLabel, rsLabel, bsLabel := getAllLabels(label)

  setsToPlotRaw := []int{}
  plotRaw(
    setsToPlotRaw,
    ras, rasLabel,
    bas, basLabel,
    rs, rsLabel,
    bs, bsLabel,
  )

  s, as := subtractBackground(ras, bas, rs, bs)

  setsToPlotSubtracted := []int{}
  plotSubtracted(setsToPlotSubtracted, s, as, rsLabel, rasLabel)

  setsToPlotSubtractedTogether := []int{}
  plotSubtractedTogether(
  setsToPlotSubtractedTogether,
  s, rsLabel,
  as, rasLabel,
  )

  subtractedGrouped := []int{}
  if len(subtractedGrouped) > 0 {
    goPlotSubGrpd(subtractedGrouped, s, as, rsLabel, rasLabel)
  }

  // Lorentz fit
  fitSets := []int{0,1,2}
  if len(fitSets) > 0 {

    // Fit parameter guesses
    amp := 12.
    wid := 0.1
    cen := 1.18

    var asAmps []float64
    var asLinewidths []float64

    fitAntiStokes := []int{}
    if len(fitAntiStokes) > 0 {

      if len(fitAntiStokes) > len(fitSets) {
        panic("fitAntiStokes must not be greater than fitSets")
      }

      // as
      fmt.Println("\nAnti-Stokes\n")
      var asFit [][]float64
      var asFits [][][]float64
      var asWidthLine [][]float64
      var asWidthLines [][][]float64
      var asfwhm []float64

      for key, set := range fitAntiStokes {

        // Test //
        // Feed previous width into next set guess -- no change
        /*if key > 0 {
          wid = asLinewidths[key - 1]/1000
          }*/

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

        // For height ratios
        asAmps = append(asAmps, amp)

        // For linewidths
        asLinewidths = append(asLinewidths, asfwhm[key])

        asFit = [][]float64{as[set][0], asyFits}
        asFits = append(asFits, asFit)
      }

      // goPlot as fits
      goPlotasFits(fitAntiStokes, as, asFits, asWidthLines, rasLabel, asfwhm, asNotes)

      // goPlot power vs width
      goPlotasPowerVsWid(fitAntiStokes, rasLabel, asNotes, asfwhm, temp)
    }

    fitStokes := []int{0,1,2}
    if len(fitStokes) > 0 {

      if len(fitStokes) > len(fitSets) {
        panic("fitStokes must not be greater than fitSets")
      }

      fmt.Println("\nStokes\n")

      var sFit [][]float64
      var sFits [][][]float64
      var sWidthLine [][]float64
      var sWidthLines [][][]float64
      var ampRatios []float64
      var sLinewidths []float64
      var sfwhm []float64

      for key, set := range fitStokes {

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

        if len(fitStokes) == len(fitAntiStokes) {
          // For height ratio
          ampRatios = append(ampRatios, amp/asAmps[key])
        }

        // For linewidth
        sLinewidths = append(sLinewidths, sfwhm[key])

        sFit = [][]float64{s[set][0], syFits}
        sFits = append(sFits, sFit)
      }
      fmt.Printf("\n")

      goPlotsFits(fitStokes, s, sFits, sWidthLines, rsLabel, sfwhm, sNotes, temp)

      goPlotsPowerVsWid(fitStokes, rsLabel, sNotes, sfwhm, temp)

      eq := true
      if len(fitAntiStokes) != len(fitStokes) {
        eq = false
      } else {
        for i, v := range fitAntiStokes {
          if v != fitStokes[i] {
            eq = false
            break
          }
        }
      }
      if eq {
        goPlotHeightRatios(fitStokes, ampRatios, asNotes, rsLabel)
        goPlotLinewidths(fitStokes, asLinewidths, sLinewidths, asNotes, sNotes, rsLabel)
      } else {
        fmt.Println("Stokes & AntiStokes sets not equal")
        fmt.Println("(Height ratio and linewidth plots not produced)\n")
      }
    }
  }
}

//--------------------------------------------------------//

func readMeta(temp bool) (string, string, []string, []string, []float64, []float64) {

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

  var date, run string
  var label, filepath []string
  var asNotes, sNotes []float64
  var dateCol, runCol, labelCol, filepathCol, notesCol int

  for col, heading := range meta[0] {
    switch heading {
    case "Date":
      dateCol = col
    case "Run":
      runCol = col
    case "Label":
      labelCol = col
    case "Filepath":
      filepathCol = col
    case "Notes":
      notesCol = col
    }
  }

  for row, value := range meta {

    if row > 0 {

      if row < 2 {
        date = value[dateCol]
        run = value[runCol]
      }
      label = append(label, value[labelCol])
      filepath = append(filepath, value[filepathCol])

      if temp {
        if strings.Contains(value[1], "as") {
          if asNote, err := strconv.ParseFloat(value[notesCol], 64); err == nil {
            asNotes = append(asNotes, asNote)
          } else {
            panic(err)
          }
        } else {
          if sNote, err := strconv.ParseFloat(value[notesCol], 64); err == nil {
            sNotes = append(sNotes, sNote)
          } else {
            panic(err)
          }
        }
      }
    }
  }

  return date, run, label, filepath, asNotes, sNotes
}

func getAllData(fileNames []string, labels []string) ([][][]float64, [][][]float64, [][][]float64, [][][]float64) {

  var bas, bs, ras, rs [][][]float64

  // Assign data by name
  for i, fileName := range fileNames {
    if strings.Contains(labels[i], "bas") {
      bas = append(bas, getData(fileName))
    } else if strings.Contains(labels[i], "bs") {
      bs = append(bs, getData(fileName))
    } else if strings.Contains(labels[i], "ras") {
      ras = append(ras, getData(fileName))
    } else if strings.Contains(labels[i], "rs") {
      rs = append(rs, getData(fileName))
    }
  }

  return ras, bas, rs, bs
}

func getData(csvName string) ([][]float64) {

  // Read
  f, err := os.Open("Data/" + csvName)
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

  for i := 1; i < len(dataStr); i++ {
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

  var rasLabel, basLabel, rsLabel, bsLabel []string

  // Assign labels by checking label
  for _, thisLabel := range label {
    if strings.Contains(thisLabel, "ras") {
      rasLabel = append(rasLabel, thisLabel)
    } else if strings.Contains(thisLabel, "bas") {
      basLabel = append(basLabel, thisLabel)
    } else if strings.Contains(thisLabel, "rs") {
      rsLabel = append(rsLabel, thisLabel)
    } else if strings.Contains(thisLabel, "bs") {
      bsLabel = append(bsLabel, thisLabel)
    }
  }

  return rasLabel, basLabel, rsLabel, bsLabel
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
  ras [][][]float64, rasLabel []string,
  bas [][][]float64, basLabel []string,
  rs [][][]float64, rsLabel []string,
  bs [][][]float64, bsLabel []string,
  ) {

  for i := 0; i < len(sets); i++ {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Raw")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    //plot.AddPointGroup(rasLabel[sets[i]], "points", ras[sets[i]])
    //plot.AddPointGroup(basLabel[0], "points", bas[0])
    plot.AddPointGroup(rsLabel[sets[i]], "points", rs[sets[i]])
    plot.AddPointGroup(bsLabel[0], "points", bs[0])
  }
}

func subtractBackground(
  ras [][][]float64,
  bas [][][]float64,
  rs [][][]float64,
  bs [][][]float64,
  ) (
  [][][]float64,
  [][][]float64,
  ) {

  var s, as [][][]float64

  for i := range rs {
    s = append(s, subtract(bs[0], rs[i]))
  }

  for i := range ras {
    as = append(as, subtract(bas[0], ras[i]))
  }

  return s, as
}

func subtract(b [][]float64, s [][]float64) ([][]float64) {

  /*var shiftUp float64 = 0

  if (s[1][0] - b[1][0] > 0) {
    shiftUp = b[1][0] - s[1][0]
  } else {
    shiftUp = s[1][600] - b[1][600]
  }*/

  for i := 0; i < len(b[0]); i++ {
    s[1][i] = s[1][i] - b[1][i] //+ shiftUp
  }

  return s
}

func plotSubtracted(
  sets []int,
  s, as [][][]float64,
  sLabel, asLabel []string,
  ) {

  for i := 0; i < len(sets); i++ {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Background Subtracted")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    plot.AddPointGroup(strings.Trim(sLabel[sets[i]], " rs") + " s", "points", s[sets[i]])
    plot.AddPointGroup(strings.Trim(asLabel[sets[i]], " pras") + " as", "points", as[sets[i]])
  }
}

func plotSubtractedTogether(
  sets []int,
  s [][][]float64, sLabel []string,
  as [][][]float64, asLabel []string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  for i := 0; i < len(sets); i++ {
    plot.AddPointGroup(strings.Trim(sLabel[sets[i]], " prs") + " s", "points", s[sets[i]])
    plot.AddPointGroup(strings.Trim(asLabel[sets[i]], " pras") + " as", "points", as[sets[i]])
  }
}

func goPlotSubGrpd(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  // as
  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Anti-Stokes"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 1
  p.X.Max = 1.3
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 36
  p.X.Tick.Label.Font.Variant = "Sans"

  /*p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 2, Label: "2"},
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
  })*/
  p.X.Padding = vg.Points(-8.5)

  p.Y.Label.Text = "Spectral Density (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 0
  p.Y.Max = 17.5
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 0, Label: "0"},
    {Value: 2.5, Label: ""},
    {Value: 5, Label: "5"},
    {Value: 7.5, Label: ""},
    {Value: 10, Label: "10"},
    {Value: 12.5, Label: ""},
    {Value: 15, Label: "15"},
    {Value: 17.5, Label: ""},
  })
  p.Y.Padding = vg.Points(-3.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)
  p.Legend.Add("Pump")

  setColors := make([]color.RGBA, 16)
  setColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
  setColors[4] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
  setColors[8] = color.RGBA{R: 122, G: 156, B: 255, A: 255}
  setColors[12] = color.RGBA{R: 255, G: 193, B: 122, A: 255}
  setColors[15] = color.RGBA{R: 27, G: 150, B: 146, A: 255}
  setColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  setColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  setColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  setColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  setColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  setColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  setColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  setColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  setColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  setColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  setColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

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

    l.GlyphStyle.Color = setColors[set]
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}
    p.Legend.Add(strings.Trim(asLabel[set], " pras"), l)
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

  savePlotAs := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/Anti-Stokes No Probe Background Subtracted"
  // Save the plot to a PNG file.
  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".png"); err != nil {
    panic(err)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".svg"); err != nil {
    panic(err)
  }

  // s
  p = plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Stokes"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 1
  p.X.Max = 1.3
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 36
  p.X.Tick.Label.Font.Variant = "Sans"

  /*p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 2, Label: "2"},
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
  })*/
  p.X.Padding = vg.Points(-8.5)

  p.Y.Label.Text = "Spectral Density (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 0
  p.Y.Max = 17.5
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 0, Label: "0"},
    {Value: 2.5, Label: ""},
    {Value: 5, Label: "5"},
    {Value: 7.5, Label: ""},
    {Value: 10, Label: "10"},
    {Value: 12.5, Label: ""},
    {Value: 15, Label: "15"},
    {Value: 17.5, Label: ""},
  })
  p.Y.Padding = vg.Points(-3.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)
  p.Legend.Add("Pump")

  for _, set := range sets {

    sPts := buildData(s[set])

    // Make a scatter plotter and set its style.
    plotSet, err := plotter.NewScatter(sPts)
    if err != nil {
      panic(err)
    }

    plotSet.GlyphStyle.Color = setColors[set]
    plotSet.GlyphStyle.Radius = vg.Points(3)
    plotSet.Shape = draw.CircleGlyph{}

    p.Add(plotSet)

    // Legend
    l, err := plotter.NewScatter(sPts)
    if err != nil {
      panic(err)
    }

    l.GlyphStyle.Color = setColors[set]
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}

    p.Legend.Add(strings.Trim(sLabel[set], " rs"), l)
  }

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

  savePlotAs = "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/Stokes No Probe Background Subtracted"
  // Save the plot to a PNG file.
  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".png"); err != nil {
    panic(err)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".svg"); err != nil {
    panic(err)
  }
}

func goPlotasFits(
  sets []int,
  as, fits, widthLines [][][]float64,
  labels []string,
  widths, notes []float64,
  ) {

    p := plot.New()
    p.BackgroundColor = color.RGBA{A:0}
    p.Title.Text = "Anti-Stokes"
    p.Title.TextStyle.Font.Typeface = "liberation"
    p.Title.TextStyle.Font.Variant = "Sans"
    p.Title.TextStyle.Font.Size = 50
    p.Title.Padding = font.Length(50)

    p.X.Label.Text = "Frequency (GHz)"
    p.X.Label.TextStyle.Font.Variant = "Sans"
    p.X.Label.TextStyle.Font.Size = 36
    p.X.Label.Padding = font.Length(20)
    p.X.LineStyle.Width = vg.Points(1.5)
    p.X.Min = 1
    p.X.Max = 1.36
    p.X.Tick.LineStyle.Width = vg.Points(1.5)
    p.X.Tick.Label.Font.Size = 36
    p.X.Tick.Label.Font.Variant = "Sans"

    /*p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: 2, Label: "2"},
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
    })*/
    p.X.Padding = vg.Points(-12.5)

    p.Y.Label.Text = "Spectral Density (nV)"
    p.Y.Label.TextStyle.Font.Variant = "Sans"
    p.Y.Label.TextStyle.Font.Size = 36
    p.Y.Label.Padding = font.Length(20)
    p.Y.LineStyle.Width = vg.Points(1.5)
    p.Y.Min = 0
    p.Y.Max = 17.5
    p.Y.Tick.LineStyle.Width = vg.Points(1.5)
    p.Y.Tick.Label.Font.Size = 36
    p.Y.Tick.Label.Font.Variant = "Sans"
    p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: 0, Label: "0"},
      {Value: 2.5, Label: ""},
      {Value: 5, Label: "5"},
      {Value: 7.5, Label: ""},
      {Value: 10, Label: "10"},
      {Value: 12.5, Label: ""},
      {Value: 15, Label: "15"},
      {Value: 17.5, Label: ""},
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

    setPtColors := make([]color.RGBA, 16)
    setPtColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
    setPtColors[4] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
    setPtColors[8] = color.RGBA{R: 122, G: 156, B: 255, A: 255}
    setPtColors[12] = color.RGBA{R: 255, G: 193, B: 122, A: 255}
    setPtColors[15] = color.RGBA{R: 27, G: 150, B: 146, A: 255}
    setPtColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
    setPtColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
    setPtColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
    setPtColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
    setPtColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
    setPtColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
    setPtColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
    setPtColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
    setPtColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
    setPtColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
    setPtColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

    setFitColors := make([]color.RGBA, 16)
    setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
    setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
    setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
    setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
    setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
    setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
    setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
    setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
    setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
    setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
    setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
    setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
    setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
    setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}


    for key, set := range sets {

      pts := buildData(as[set])
      fit := buildData(fits[key])
      wid := buildData(widthLines[key])

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
      power := strings.Trim(labels[set], " pras")
      temp := strconv.FormatFloat(notes[set], 'f', -1, 64)
      p.Legend.Add(power + " @" + temp + "K", l)
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

func goPlotasPowerVsWid(
  sets []int,
  labels []string,
  notes, widths []float64,
  temp bool,
  ) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Anti-Stokes Pump Power vs Widths of Fits"
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
  /*p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 90, Label: "90"},
    {Value: 95, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 105, Label: ""},
    {Value: 110, Label: "110"},
    {Value: 115, Label: ""},
    {Value: 120, Label: "120"},
    {Value: 125, Label: ""},
    {Value: 130, Label: "130"},
  })*/
  p.Y.Padding = vg.Points(1)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setFitColors := make([]color.RGBA, 16)
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
  setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
  setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

  for key, set := range sets {

    pts := make(plotter.XYs, 1)

    power := strings.Trim(labels[set], " mW pras")
    if pwr, err := strconv.ParseFloat(power, 64); err == nil {
      pts[0].X = pwr
    } else {
      panic(err)
    }
    pts[0].Y = widths[key]

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
    if temp {
      temperature := strconv.FormatFloat(notes[set], 'f', -1, 64)
      p.Legend.Add(power + " mW @" + temperature + "K", plotPts)
    } else {
      p.Legend.Add(power + " mW")
    }
  }

  // Save plot
  name := "as Pow vs Wid"
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

func goPlotsFits(
  sets []int,
  s, fits, widthLines [][][]float64,
  labels []string,
  widths, notes []float64,
  temp bool,
  ) {

    p := plot.New()
    p.BackgroundColor = color.RGBA{A:0}
    p.Title.Text = "Stokes"
    p.Title.TextStyle.Font.Typeface = "liberation"
    p.Title.TextStyle.Font.Variant = "Sans"
    p.Title.TextStyle.Font.Size = 50
    p.Title.Padding = font.Length(50)

    p.X.Label.Text = "Frequency (GHz)"
    p.X.Label.TextStyle.Font.Variant = "Sans"
    p.X.Label.TextStyle.Font.Size = 36
    p.X.Label.Padding = font.Length(20)
    p.X.LineStyle.Width = vg.Points(1.5)
    p.X.Min = 1
    p.X.Max = 1.36
    p.X.Tick.LineStyle.Width = vg.Points(1.5)
    p.X.Tick.Label.Font.Size = 36
    p.X.Tick.Label.Font.Variant = "Sans"

    p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: 1, Label: "1"},
      {Value: 1.05, Label: ""},
      {Value: 1.1, Label: "1.1"},
      {Value: 1.15, Label: ""},
      {Value: 1.2, Label: "1.2"},
      {Value: 1.25, Label: ""},
      {Value: 1.3, Label: "1.3"},
      {Value: 1.35, Label: ""},
      {Value: 1.4, Label: "1.4"},
    })
    p.X.Padding = vg.Points(-12.5)

    p.Y.Label.Text = "Spectral Density (nV)"
    p.Y.Label.TextStyle.Font.Variant = "Sans"
    p.Y.Label.TextStyle.Font.Size = 36
    p.Y.Label.Padding = font.Length(20)
    p.Y.LineStyle.Width = vg.Points(1.5)
    p.Y.Min = 0
    p.Y.Max = 17.5
    p.Y.Tick.LineStyle.Width = vg.Points(1.5)
    p.Y.Tick.Label.Font.Size = 36
    p.Y.Tick.Label.Font.Variant = "Sans"
    p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
      {Value: 0, Label: "0"},
      {Value: 2.5, Label: ""},
      {Value: 5, Label: "5"},
      {Value: 7.5, Label: ""},
      {Value: 10, Label: "10"},
      {Value: 12.5, Label: ""},
      {Value: 15, Label: "15"},
      {Value: 17.5, Label: ""},
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

    setPtColors := make([]color.RGBA, 16)
    setPtColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
    setPtColors[4] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
    setPtColors[8] = color.RGBA{R: 122, G: 156, B: 255, A: 255}
    setPtColors[12] = color.RGBA{R: 255, G: 193, B: 122, A: 255}
    setPtColors[15] = color.RGBA{R: 27, G: 150, B: 146, A: 255}
    setPtColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
    setPtColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
    setPtColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
    setPtColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
    setPtColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
    setPtColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
    setPtColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
    setPtColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
    setPtColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
    setPtColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
    setPtColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

    setFitColors := make([]color.RGBA, 16)
    setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
    setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
    setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
    setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
    setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
    setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
    setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
    setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
    setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
    setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
    setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
    setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
    setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
    setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}


    for key, set := range sets {

      pts := buildData(s[set])
      fit := buildData(fits[key])
      wid := buildData(widthLines[key])

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
      power := strings.Trim(labels[set], " prs")
      if temp {
        temperature := strconv.FormatFloat(notes[set], 'f', -1, 64)
        p.Legend.Add(power + " @" + temperature + "K", l)
      } else {
        p.Legend.Add(power)
      }
    }

    // Save plot
    name := "Stokes w Fits"
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

func goPlotsPowerVsWid(
  sets []int,
  labels []string,
  notes, widths []float64,
  temp bool,
  ) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Stokes Pump Power vs Widths of Fits"
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
  p.Y.Min = 0
  p.Y.Max = 500
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 0, Label: "0"},
    {Value: 50, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 150, Label: ""},
    {Value: 200, Label: "200"},
    {Value: 250, Label: ""},
    {Value: 300, Label: "300"},
    {Value: 350, Label: ""},
    {Value: 400, Label: "400"},
    {Value: 450, Label: ""},
    {Value: 500, Label: "500"},
  })
  p.Y.Padding = vg.Points(1)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setFitColors := make([]color.RGBA, 16)
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
  setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
  setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

  for key, set := range sets {

    pts := make(plotter.XYs, 1)

    power := strings.Trim(labels[set], " mW prs")
    if pwr, err := strconv.ParseFloat(power, 64); err == nil {
      pts[0].X = pwr
    } else {
      panic(err)
    }
    pts[0].Y = widths[key]

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
    v[0].Y = 0
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
    if temp {
      temperature := strconv.FormatFloat(notes[set], 'f', -1, 64)
      p.Legend.Add(power + " mW @" + temperature + "K", plotPts)
    } else {
      p.Legend.Add(power + " mW")
    }
  }

  // Save plot
  name := "s Pow vs Wid"
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

func goPlotHeightRatios(sets []int, heightRatios, powers []float64, labels []string) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Height Ratios vs Power"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Pump Power (mW)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 25
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

  p.Y.Label.Text = "Stokes/Anti-Stokes Heights"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 1
  p.Y.Max = 1.3
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  /*p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 90, Label: "90"},
    {Value: 95, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 105, Label: ""},
    {Value: 110, Label: "110"},
    {Value: 115, Label: ""},
    {Value: 120, Label: "120"},
    {Value: 125, Label: ""},
    {Value: 130, Label: "130"},
  })*/
  p.Y.Padding = vg.Points(-.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setFitColors := make([]color.RGBA, 16)
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
  setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
  setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

  // Linear fit line

  // Fit parameter guesses
  m := 5.
  b := .0125

  f := func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for key, _ := range sets {

      x = powers[key]
      y := heightRatios[key]

      dst[key] = m * x + b - y
    }
  }

  jacobian := lm.NumJac{Func: f}

  // Solve for fit
  toBeSolved := lm.LMProblem{
    Dim:        2,
    Size:       len(sets),
    Func:       f,
    Jac:        jacobian.Jac,
    InitParams: []float64{m, b},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }

  results, _ := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  m, b = results.X[0], results.X[1]

  var yFit []float64
  var xFit []float64

  // Create function according to solved fit parameters
  for key, _ := range sets {
    var x float64

    x = powers[key]

    xFit = append(xFit, x)
    yFit = append(yFit, m * x + b)
  }

  fit := buildData([][]float64{xFit, yFit})

  // Plot fit
  plotFit, err := plotter.NewLine(fit)
  if err != nil {
    panic(err)
  }

  p.Add(plotFit)

  plotFit.LineStyle.Color = color.RGBA{R: 127, G: 127, B: 127, A: 255}
  plotFit.LineStyle.Width = vg.Points(3)

  for key, _ := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = powers[key]
    pts[0].Y = heightRatios[key]

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      panic(err)
    }

    plotPts.GlyphStyle.Color = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    plotPts.GlyphStyle.Radius = vg.Points(6)
    plotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(plotPts)
    //p.Legend.Add(strings.Trim(labels[set], " prs"), plotPts)
  }

  // Save plot
  name := "height ratios"
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

func goPlotLinewidths(sets []int, asLinewidths, sLinewidths, asPowers, sPowers []float64, labels []string) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = "Linewidths vs Power"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Pump Power (mW)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 25
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
  p.Y.Min = 50
  p.Y.Max = 125
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 50, Label: "50"},
    {Value: 62.5, Label: ""},
    {Value: 75, Label: "75"},
    {Value: 87.5, Label: ""},
    {Value: 100, Label: "100"},
    {Value: 112.5, Label: ""},
    {Value: 125, Label: "125"},
  })
  p.Y.Padding = vg.Points(-.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setFitColors := make([]color.RGBA, 16)
  setFitColors[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
  setFitColors[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  setFitColors[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  setFitColors[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
  setFitColors[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
  setFitColors[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  setFitColors[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  setFitColors[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  setFitColors[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  setFitColors[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  setFitColors[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  setFitColors[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  setFitColors[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  setFitColors[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  setFitColors[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  setFitColors[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

  // as linear fit

  // Fit parameter guesses
  m := 1.
  b := 100.

  f := func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for key, _ := range sets {

      x = asPowers[key]
      y := asLinewidths[key]

      dst[key] = m * x + b - y
    }
  }

  jacobian := lm.NumJac{Func: f}

  // Solve for fit
  toBeSolved := lm.LMProblem{
    Dim:        2,
    Size:       len(sets),
    Func:       f,
    Jac:        jacobian.Jac,
    InitParams: []float64{m, b},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }

  results, _ := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  m, b = results.X[0], results.X[1]

  var asyFit []float64
  var asxFit []float64

  // Create function according to solved fit parameters
  for key, _ := range sets {
    var x float64

    x = asPowers[key]

    asxFit = append(asxFit, x)
    asyFit = append(asyFit, m * x + b)
  }

  asfit := buildData([][]float64{asxFit, asyFit})

  // Plot as fit
  asPlotFit, err := plotter.NewLine(asfit)
  if err != nil {
    panic(err)
  }

  asPlotFit.LineStyle.Color = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  asPlotFit.LineStyle.Width = vg.Points(3)

  // s linear fit

  f = func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for key, _ := range sets {

      x = sPowers[key]
      y := sLinewidths[key]

      dst[key] = m * x + b - y
    }
  }

  jacobian = lm.NumJac{Func: f}

  // Solve for fit
  toBeSolved = lm.LMProblem{
    Dim:        2,
    Size:       len(sets),
    Func:       f,
    Jac:        jacobian.Jac,
    InitParams: []float64{m, b},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }

  results, _ = lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  m, b = results.X[0], results.X[1]

  var syFit []float64
  var sxFit []float64

  // Create function according to solved fit parameters
  for key, _ := range sets {
    var x float64

    x = sPowers[key]

    sxFit = append(sxFit, x)
    syFit = append(syFit, m * x + b)
  }

  sfit := buildData([][]float64{sxFit, syFit})

  // Plot as fit
  sPlotFit, err := plotter.NewLine(sfit)
  if err != nil {
    panic(err)
  }

  sPlotFit.LineStyle.Color = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  sPlotFit.LineStyle.Width = vg.Points(3)

  p.Add(asPlotFit, sPlotFit)
  p.Legend.Add("Anti-Stokes", asPlotFit)
  p.Legend.Add("Stokes", sPlotFit)

  // as points
  for key, _ := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = asPowers[key]
    pts[0].Y = asLinewidths[key]

    // Plot points
    asPlotPts, err := plotter.NewScatter(pts)
    if err != nil {
      panic(err)
    }

    asPlotPts.GlyphStyle.Color = color.RGBA{R: 99, G: 124, B: 198, A: 255}
    asPlotPts.GlyphStyle.Radius = vg.Points(6)
    asPlotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(asPlotPts)
  }

  // s points
  for key, _ := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = sPowers[key]
    pts[0].Y = sLinewidths[key]

    // Plot points
    sPlotPts, err := plotter.NewScatter(pts)
    if err != nil {
      panic(err)
    }

    sPlotPts.GlyphStyle.Color = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    sPlotPts.GlyphStyle.Radius = vg.Points(6)
    sPlotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(sPlotPts)
  }

  // Save plot
  name := "linewidths"
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
