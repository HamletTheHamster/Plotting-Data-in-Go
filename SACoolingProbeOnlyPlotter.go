package main

import (
  "image/color"
  "github.com/Arafatk/glot"
  //"github.com/maorshutman/lm"
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
  plotSubtracted(
    setsToPlotSubtracted,
    s, rsLabel,
    as, rasLabel,
  )

  setsToPlotSubtractedTogether := []int{}
  plotSubtractedTogether(
  setsToPlotSubtractedTogether,
  s, rsLabel,
  as, rasLabel,
  )

  subtractedGrouped := []int{0,1,2,3,4}
  if len(subtractedGrouped) > 0 {
    goPlotSubGrpd(subtractedGrouped, s, as, rsLabel, rasLabel)
  }

  // Lorentz Fit

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
  f, err := os.Open(csvName)
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
    var convSignal []float64

    for _, sigElemToConvert := range signal {
      convSignal = append(convSignal, math.Pow(10, 6)*math.Pow(10, sigElemToConvert/10.))
    }

    return [][]float64{frequency, convSignal}
  } else if dataStr[1][3] == "  uV" {

    return [][]float64{frequency, signal}
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

    plot.AddPointGroup(rasLabel[sets[i]], "points", ras[sets[i]])
    plot.AddPointGroup(basLabel[0], "points", bas[0])
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

  for i := 0; i < len(ras); i++ {
    s = append(s, subtract(bs[0], rs[i]))
    as = append(as, subtract(bas[0], ras[i]))
  }

  return s, as
}

func subtract(b [][]float64, s [][]float64) ([][]float64) {

  var shiftUp float64 = 0

  if (s[1][0] - b[1][0] > 0) {
    shiftUp = b[1][0] - s[1][0]
  }

  for i := 0; i < len(b[0]); i++ {
    s[1][i] = s[1][i] - b[1][i] + shiftUp
  }

  return s
}

func plotSubtracted(
  sets []int,
  s [][][]float64, sLabel []string,
  as [][][]float64, asLabel []string,
  ) {

  for i := 0; i < len(sets); i++ {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Background Subtracted")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    plot.AddPointGroup(strings.Trim(sLabel[sets[i]], " prs") + " s", "points", s[sets[i]])
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
  p.Title.Text = "Anti-Stokes No Probe"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 34
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 24
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 2.1
  p.X.Max = 2.4
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 24
  p.X.Tick.Label.Font.Variant = "Sans"

  p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 2.15, Label: ""},
    {Value: 2.2, Label: "2.2"},
    {Value: 2.25, Label: ""},
    {Value: 2.3, Label: "2.3"},
    {Value: 2.35, Label: ""},
    {Value: 2.4, Label: "2.4"},
  })
  p.X.Padding = vg.Points(-8.5)

  p.Y.Label.Text = "Signal (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 24
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 0
  p.Y.Max = 0.75
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 24
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: .05, Label: ""},
    {Value: .1, Label: ".1"},
    {Value: .15, Label: ""},
    {Value: .2, Label: ".2"},
    {Value: .25, Label: ""},
    {Value: .3, Label: ".3"},
    {Value: .35, Label: ""},
    {Value: .4, Label: ".4"},
    {Value: .45, Label: ""},
    {Value: .5, Label: ".5"},
    {Value: .55, Label: ""},
    {Value: .6, Label: ".6"},
    {Value: .65, Label: ""},
    {Value: .7, Label: ".7"},
    {Value: .75, Label: ""},
  })
  p.Y.Padding = vg.Points(-3.75)

  p.Legend.TextStyle.Font.Size = 24
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setColors := make([]color.RGBA, 5)
  setColors[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
  setColors[1] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
  setColors[2] = color.RGBA{R: 122, G: 156, B: 255, A: 255}
  setColors[3] = color.RGBA{R: 255, G: 243, B: 117, A: 255}
  setColors[4] = color.RGBA{R: 188, G: 117, B: 255, A: 255}


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
    p.Legend.Add(strings.Trim(asLabel[set], " pras"), plotSet)
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
  p.Title.Text = "Stokes No Probe"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 34
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 24
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = 2.1
  p.X.Max = 2.4
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 24
  p.X.Tick.Label.Font.Variant = "Sans"

  p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: 2.15, Label: ""},
    {Value: 2.2, Label: "2.2"},
    {Value: 2.25, Label: ""},
    {Value: 2.3, Label: "2.3"},
    {Value: 2.35, Label: ""},
    {Value: 2.4, Label: "2.4"},
  })
  p.X.Padding = vg.Points(-8.5)

  p.Y.Label.Text = "Signal (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 24
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = 0
  p.Y.Max = 0.75
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 24
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
    {Value: .05, Label: ""},
    {Value: .1, Label: ".1"},
    {Value: .15, Label: ""},
    {Value: .2, Label: ".2"},
    {Value: .25, Label: ""},
    {Value: .3, Label: ".3"},
    {Value: .35, Label: ""},
    {Value: .4, Label: ".4"},
    {Value: .45, Label: ""},
    {Value: .5, Label: ".5"},
    {Value: .55, Label: ""},
    {Value: .6, Label: ".6"},
    {Value: .65, Label: ""},
    {Value: .7, Label: ".7"},
    {Value: .75, Label: ""},
  })
  p.Y.Padding = vg.Points(-3.75)

  p.Legend.TextStyle.Font.Size = 24
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

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
    p.Legend.Add(strings.Trim(sLabel[set], " rs"), plotSet)
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

func normalizeFit(fit []float64) ([]float64) {

  var shift float64 = (fit[0] + fit[599])/2

  for i := 0; i < 600; i++ {
    fit[i] = fit[i] - shift
  }
  return fit
}
