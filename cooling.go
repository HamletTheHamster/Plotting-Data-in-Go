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

  lock, temp, lcof, uhna3 := flags()

  date, run, label, asPowers, sPowers, file, asNotes, sNotes := readMeta(lock, temp)
  ras, bas, rs, bs := getAllData(lock, file, label)

  log := header(lock, temp, lcof, date, run)

  asLabel, basLabel, sLabel, bsLabel := getAllLabels(label)

  setsToPlotRaw := []int{}
  plotRaw(
    setsToPlotRaw,
    bas, ras, bs, rs,
    basLabel, asLabel, bsLabel, sLabel,
  )

  s, as := subtractBackground(ras, bas, rs, bs)

  setsToPlotSubtracted := []int{}
  plotSubtracted(
    setsToPlotSubtracted,
    s, as,
    sLabel, asLabel,
  )

  setsToPlotSubtractedTogether := []int{}
  plotSubtractedTogether(
    setsToPlotSubtractedTogether,
    as, s,
    asLabel, sLabel,
  )

  subtractedGrouped := []int{}
  if len(subtractedGrouped) > 0 {
    goPlotSubGrpd(subtractedGrouped, s, as, sLabel, asLabel)
  }

  fitSets := true
  if fitSets {

    var amp, wid, cen float64
    var sample string

    if lcof {
      sample = "Liquid-Core"
      amp = 5.
      wid = 0.1
      cen = 2.25
    } else if uhna3 {
      sample = "UHNA3"
      amp = 12.
      wid = 0.1
      cen = 1.18
    } else {
      sample = "[Unspecified] Sample"
      amp = 5.
      wid = 0.1
      cen = 2.25
    }

    var asAmps, asLinewidths []float64

    // 0,4,8,12,15: presentation
    // 2,3,4,5,6,7,8,9: <
    fitAntiStokes := []int{0,2,7,11,15}
    if len(fitAntiStokes) > 0 {

      // as
      header := fmt.Sprintf("\nAnti-Stokes\nSet  \t Power \t\t Width \t\t Peak \t\t Center\n")
      fmt.Printf(header)
      log = append(log, header)

      var asFit [][]float64
      var asFits [][][]float64
      var asWidthLine [][]float64
      var asWidthLines [][][]float64
      var asfwhm []float64

      for i, set := range fitAntiStokes {

        f := func(dst, guess []float64) {

          amp, wid, cen := guess[0], guess[1], guess[2]

          for i := range as[set][0] {
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

        str := fmt.Sprintf("%d \t %.2f mW \t %.2f MHz \t %.6f nV \t %.4f GHz\n", set, asPowers[set], wid*2000, amp, cen)
        fmt.Printf(str)
        log = append(log, str)

        var asyFits []float64

        // Create function according to solved fit parameters
        for i := range as[set][0] {
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
        asLinewidths = append(asLinewidths, asfwhm[i])

        asFit = [][]float64{as[set][0], asyFits}
        asFits = append(asFits, asFit)
      }

      // goPlot as fits
      goPlotasFits(fitAntiStokes, as, asFits, asWidthLines, asLabel, asfwhm, asNotes, temp, sample)

      // goPlot power vs width
      goPlotasPowerVsWid(fitAntiStokes, asLabel, asNotes, asfwhm, temp, sample)
    }

    fitStokes := []int{0,2,7,11,15}
    if len(fitStokes) > 0 {

      header := "\nStokes\nSet \t Power \t\t Width \t\t Peak \t\t Center \n"
      fmt.Printf(header)
      log = append(log, header)

      var sFit [][]float64
      var sFits [][][]float64
      var sWidthLine [][]float64
      var sWidthLines [][][]float64
      var ampRatios []float64
      var sLinewidths []float64
      var sfwhm []float64

      for i, set := range fitStokes {

        f := func(dst, guess []float64) {

          amp, wid, cen := guess[0], guess[1], guess[2]

          for i := range s[set][0] {
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

        str := fmt.Sprintf("%d \t %.2f mW \t %.2f MHz \t %.6f nV \t %.4f GHz\n", set, sPowers[set], wid*2000, amp, cen)
        fmt.Printf(str)
        log = append(log, str)

        var syFits []float64

        // Create function according to solved fit parameters
        for i := range s[set][0] {
          // (amp*wid^2/((x-cen)^2+wid^2))
          x := s[set][0][i]
          syFits = append(syFits, results.X[0] * math.Pow(results.X[1], 2) / (math.Pow(x - results.X[2], 2) + math.Pow(results.X[1], 2)))
        }

        // Width lines
        sWidthLine = [][]float64{{cen - wid, cen + wid},{amp/2, amp/2}}
        sWidthLines = append(sWidthLines, sWidthLine)

        if len(fitStokes) == len(fitAntiStokes) {
          // For height ratio
          ampRatios = append(ampRatios, amp/asAmps[i])
        }

        // For linewidth
        sLinewidths = append(sLinewidths, sfwhm[i])

        sFit = [][]float64{s[set][0], syFits}
        sFits = append(sFits, sFit)
      }
      fmt.Printf("\n")
      log = append(log, "\n")

      goPlotsFits(fitStokes, s, sFits, sWidthLines, sLabel, sfwhm, sNotes, temp, sample)

      goPlotsPowerVsWid(fitStokes, sLabel, sNotes, sfwhm, temp, sample)

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
        var powers []float64
        for i, v := range asPowers {
          powers = append(powers, (v + sPowers[i])/2)
        }
        goPlotHeightRatios(fitStokes, ampRatios, powers, sLabel, sample)
        goPlotLinewidths(fitStokes, asLinewidths, sLinewidths, asPowers, sPowers, sLabel, sample)
      } else {
        str := fmt.Sprintf("Stokes & AntiStokes sets not equal\n" +
          "(Height ratio and linewidth plots not produced)\n")
        fmt.Printf(str)
        log = append(log, str)
      }
    }

    // binMHz := 50.
    // as, s = bin(as,s, binMHz)
  }

  writeLog(log)
}

//----------------------------------------------------------------------------//

func flags() (
  bool, bool, bool, bool,
) {

  var lock, temp, lcof, uhna3 bool

  flag.BoolVar(&lock, "l", false, "lock-in data")
  flag.BoolVar(&temp, "t", false, "contains temperature data")
  flag.BoolVar(&lcof, "o", false, "liquid-core optical fiber sample")
  flag.BoolVar(&uhna3, "3", false, "UHNA3 fiber sample")
  flag.Parse()

  if lcof && uhna3 {
    fmt.Println("flag.Parse(): sample flagged as both UHNA3 and liquid-core.")
    os.Exit(1)
  }

  return lock, temp, lcof, uhna3
}

func readMeta(
  lock, temp bool,
) (
  string, string, []string, []float64, []float64, []string, []float64, []float64,
) {

  // Read
  metaFile, err := os.Open("Data/meta.csv")
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  reader := csv.NewReader(metaFile)
  meta, err := reader.ReadAll()
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  var date, run string
  var label, filepath []string
  var asPowers, sPowers, asNotes, sNotes []float64
  var dateCol, runCol, labelCol, filepathCol, sigFileCol, freqFileCol, notesCol int

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
    case "Signal":
      sigFileCol = col
    case "Frequency":
      freqFileCol = col
    case "Notes":
      notesCol = col
    }
  }

  for row, v := range meta {

    if row > 0 {

      if row < 2 {
        date = v[dateCol]
        run = v[runCol]
      }

      label = append(label, v[labelCol])

      if strings.Contains(v[labelCol], "ras") {
        if v, err := strconv.ParseFloat(strings.Split(v[labelCol], " ")[0], 64); err == nil {
          asPowers = append(asPowers, v)
        } else {
          fmt.Println(err)
          os.Exit(1)
        }
      } else if strings.Contains(v[labelCol], "rs"){
        if v, err := strconv.ParseFloat(strings.Split(v[labelCol], " ")[0], 64); err == nil {
          sPowers = append(sPowers, v)
        } else {
          fmt.Println(err)
          os.Exit(1)
        }
      }

      if lock {
        filepath = append(filepath, v[sigFileCol])
        filepath = append(filepath, v[freqFileCol])
      } else {
        filepath = append(filepath, v[filepathCol])
      }

      if temp {
        if strings.Contains(v[labelCol], "as") {
          if asNote, err := strconv.ParseFloat(v[notesCol], 64); err == nil {
            asNotes = append(asNotes, asNote)
          } else {
            fmt.Println(err)
            os.Exit(1)
          }
        } else {
          if sNote, err := strconv.ParseFloat(v[notesCol], 64); err == nil {
            sNotes = append(sNotes, sNote)
          } else {
            fmt.Println(err)
            os.Exit(1)
          }
        }
      }
    }
  }

  return date, run, label, asPowers, sPowers, filepath, asNotes, sNotes
}

func getAllData(
  lock bool,
  fileNames, labels []string,
) (
  [][][]float64, [][][]float64, [][][]float64, [][][]float64,
) {

  var bas, bs, ras, rs [][][]float64
  sig := false

  if lock {

    // Assign data by name
    for i := 0; i < len(fileNames)/2; i++ {

      if i == 0 || i%2 != 0 {
        sig = true
      }

      if strings.Contains(labels[i], "bas") {
        bas = append(bas, getData(lock, sig, fileNames[i]))
      } else if strings.Contains(labels[i], "bs") {
        bs = append(bs, getData(lock, sig, fileNames[i]))
      } else if strings.Contains(labels[i], "ras") {
        ras = append(ras, getData(lock, sig, fileNames[i]))
      } else if strings.Contains(labels[i], "rs") {
        rs = append(rs, getData(lock, sig, fileNames[i]))
      }
    }
  } else {

    // Assign data by name
    for i, fileName := range fileNames {

      if i == 0 || i%2 != 0 {
        sig = true
      }

      if strings.Contains(labels[i], "bas") {
        bas = append(bas, getData(lock, sig, fileName))
      } else if strings.Contains(labels[i], "bs") {
        bs = append(bs, getData(lock, sig, fileName))
      } else if strings.Contains(labels[i], "ras") {
        ras = append(ras, getData(lock, sig, fileName))
      } else if strings.Contains(labels[i], "rs") {
        rs = append(rs, getData(lock, sig, fileName))
      }
    }
  }

  return ras, bas, rs, bs
}

func getData(
  lock, sig bool,
  csvName string,
) (
  [][]float64,
) {

  // Read
  f, err := os.Open("Data/" + csvName)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }
  defer f.Close()
  dataStr, err := readCSV(f)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  // Separate, Strip, & Transpose
  var frequencyStrT, signalStrT []string

  if lock {
    if sig {
      for i := range dataStr {
        signalStrT = append(signalStrT, dataStr[i][0])
      }
    } else {
      for i := range dataStr {
        frequencyStrT = append(frequencyStrT, dataStr[i][0])
      }
    }
  } else {
    for i := 1; i < len(dataStr); i++ {
      frequencyStrT = append(frequencyStrT, strings.ReplaceAll(dataStr[i][0]," ",""))
      signalStrT = append(signalStrT, strings.ReplaceAll(dataStr[i][2]," ",""))
    }
  }

  // Convert to float
  var frequency, signal []float64

  for _, freqElem := range frequencyStrT {
    if freqValue, err := strconv.ParseFloat(freqElem, 64); err != nil {
      fmt.Println(err)
      os.Exit(1)
    } else {
      frequency = append(frequency, freqValue/1e9)
    }
  }

  for _, sigElem := range signalStrT {
    if sigValue, err := strconv.ParseFloat(sigElem, 64); err != nil {
      fmt.Println(err)
      os.Exit(1)
    } else {
      signal = append(signal, sigValue)
    }
  }

  if !lock {
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
  } else {
    return [][]float64{frequency, signal}
  }
}

func readCSV(
  rs io.ReadSeeker,
) (
  [][]string, error,
) {
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

func header(
  lock, temp, lcof bool,
  date, run string,
) (
  []string,
) {

  log := []string{}
  log = append(log, date)
  if run != "" {
    log = append(log, date + " Run " + run + "\n")
  }

  fmt.Printf(log[0])

  if temp {
    log = append(log, "\n*Temperature-dependent data*\n")
    fmt.Printf("\n*Temperature-dependent data*\n")
  }
  if lcof {
    log = append(log, "\n*Liquid-core optical fiber sample*\n")
    fmt.Printf("\n*Liquid-core optical fiber sample*\n")
  }
  if lock {
    log = append(log, "\n*Data gathered from Lock-in*\n")
    fmt.Printf("\n*Data gathered from Lock-in*\n")
  } else {
    log = append(log, "\n*Data gathered from Spectrum Analyzer*\n")
    fmt.Printf("\n*Data gathered from Spectrum Analyzer*\n")
  }

  return log
}

func getAllLabels(
  label []string,
) (
  []string, []string, []string, []string,
) {

  var rasLabel, basLabel, rsLabel, bsLabel []string

  // Assign labels by verifying label
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

func buildData(
  data [][]float64,
) (
  plotter.XYs,
) {

  xy := make(plotter.XYs, len(data[0]))

  for i := range xy {
    xy[i].X = data[0][i]
    xy[i].Y = data[1][i]
  }

  return xy
}

func plotRaw(
  sets []int,
  bas, ras, bs, rs [][][]float64,
  basLabel, rasLabel, bsLabel, rsLabel []string,
) {

  for i := range sets {
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

func axes(
  plot, sample string,
) (
  []float64, []float64, []float64, []float64, []string, []string, error,
) {

  switch plot {
  case "fits":
    switch sample {
    case "Liquid-Core":
      xrange := []float64{2, 2.5}
      yrange := []float64{0, 17.5}
      xtick := []float64{2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5}
      ytick := []float64{0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5}
      xtickLabel := []string{"2", "", "2.1", "", "2.2", "", "2.3", "", "2.4", "", "2.5"}
      ytickLabel := []string{"0", "", "5", "", "10", "", "15", ""}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "UHNA3":
      xrange := []float64{1, 1.36}
      yrange := []float64{0, 17.5}
      xtick := []float64{1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}
      ytick := []float64{0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5}
      xtickLabel := []string{"1", "", "1.1", "", "1.2", "", "1.3", "", "1.4"}
      ytickLabel := []string{"0", "", "5", "", "10", "", "15", ""}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "[Unspecified Sample]":
      xrange := []float64{2, 2.5}
      yrange := []float64{0, 17.5}
      xtick := []float64{2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5}
      ytick := []float64{0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5}
      xtickLabel := []string{"2", "", "2.1", "", "2.2", "", "2.3", "", "2.4", "", "2.5"}
      ytickLabel := []string{"0", "", "5", "", "10", "", "15", ""}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    }
  case "pow vs wid":
    switch sample {
    case "Liquid-Core":
      xrange := []float64{0, 200}
      yrange := []float64{90, 130}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{90, 95, 100, 105, 110, 115, 120, 125, 130}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"90", "", "100", "", "110", "", "120", "", "130"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "UHNA3":
      xrange := []float64{0, 200}
      yrange := []float64{90, 130}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{90, 95, 100, 105, 110, 115, 120, 125, 130}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"90", "", "100", "", "110", "", "120", "", "130"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "[Unspecified Sample]":
      xrange := []float64{0, 200}
      yrange := []float64{90, 130}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{90, 95, 100, 105, 110, 115, 120, 125, 130}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"90", "", "100", "", "110", "", "120", "", "130"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    }
  case "height ratios":
    switch sample {
    case "Liquid-Core":
      xrange := []float64{0, 200}
      yrange := []float64{1, 3.5}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{1, 1.5, 2, 2.5, 3, 3.5}
      xtickLabel := []string{"", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"1", "", "2", "", "3", ""}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "UHNA3":
      xrange := []float64{0, 200}
      yrange := []float64{90, 130}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{90, 95, 100, 105, 110, 115, 120, 125, 130}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"90", "", "100", "", "110", "", "120", "", "130"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "[Unspecified Sample]":
      xrange := []float64{0, 200}
      yrange := []float64{90, 130}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{90, 95, 100, 105, 110, 115, 120, 125, 130}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"90", "", "100", "", "110", "", "120", "", "130"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    }
  case "linewidths":
    switch sample {
    case "Liquid-Core":
      xrange := []float64{75, 175}
      yrange := []float64{62.5, 150}
      xtick := []float64{75, 100, 125, 150, 175}
      ytick := []float64{62.5, 75, 87.5, 100, 112.5, 125, 137.5, 150}
      xtickLabel := []string{"75", "", "100", "", "150", ""}
      ytickLabel := []string{"", "75", "", "100", "", "125", "", "150"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "UHNA3":
      xrange := []float64{25, 200}
      yrange := []float64{50, 125}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{50, 62.5, 75, 87.5, 100, 112.5, 125}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"50", "", "75", "", "100", "", "125"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    case "[Unspecified Sample]":
      xrange := []float64{25, 200}
      yrange := []float64{50, 125}
      xtick := []float64{0, 25, 50, 75, 100, 125, 150, 175, 200}
      ytick := []float64{50, 62.5, 75, 87.5, 100, 112.5, 125}
      xtickLabel := []string{"0", "", "50", "", "100", "", "150", "", "200"}
      ytickLabel := []string{"50", "", "75", "", "100", "", "125"}

      return xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, nil
    }
  }

  return []float64{}, []float64{}, []float64{}, []float64{}, []string{},
    []string{}, fmt.Errorf("func axes: no predefined axis for '%s'", plot)
}

func subtractBackground(
  ras, bas, rs, bs [][][]float64,
) (
  [][][]float64, [][][]float64,
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

func subtract(
  b, s [][]float64,
) (
  [][]float64,
) {

  /*var shiftUp float64 = 0

  if (s[1][0] - b[1][0] > 0) {
    shiftUp = b[1][0] - s[1][0]
  } else {
    shiftUp = s[1][600] - b[1][600]
  }*/

  for i := range b[0] {
    s[1][i] = s[1][i] - b[1][i] //+ shiftUp
  }

  return s
}

func plotSubtracted(
  sets []int,
  s, as [][][]float64,
  sLabel, asLabel []string,
) {

  for _, set := range sets {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Background Subtracted")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    plot.AddPointGroup(strings.Trim(sLabel[sets[set]], " rs") + " s", "points", s[sets[set]])
    plot.AddPointGroup(strings.Trim(asLabel[sets[set]], " ras") + " as", "points", as[sets[set]])
  }
}

func plotSubtractedTogether(
  sets []int,
  as, s [][][]float64,
  asLabel, sLabel []string,
) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  for _, set := range sets {
    plot.AddPointGroup(strings.Trim(sLabel[sets[set]], " rs") + " s", "points", s[sets[set]])
    //plot.AddPointGroup(strings.Trim(asLabel[sets[set]], " ras") + " as", "points", as[sets[set]])
  }
}

func goPlotSubGrpd(
  sets []int,
  s, as [][][]float64,
  sLabel, asLabel []string,
) {

  // Anti-Stokes
  title := "Anti-Stokes"
  xlabel := "Frequency (GHz)"
  ylabel := "Spectral Density (nV)"
  legend := "Pump"
  xrange := []float64{1, 1.36}
  yrange := []float64{0, 17.5}
  xtick := []float64{1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}
  ytick := []float64{0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5}
  xtickLabels := []string{"1", "", "1.1", "", "1.2", "", "1.3", "", "1.4"}
  ytickLabels := []string{"0", "", "5", "", "10", "", "15", ""}

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabels, ytickLabels,
  )

  for _, set := range sets {

    asPts := buildData(as[set])

    // Make a scatter plotter and set its style.
    plotSet, err := plotter.NewScatter(asPts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotSet.GlyphStyle.Color = palette(set, false)
    plotSet.GlyphStyle.Radius = vg.Points(3)
    plotSet.Shape = draw.CircleGlyph{}

    p.Add(plotSet)

    // Legend
    l, err := plotter.NewScatter(asPts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    l.GlyphStyle.Color = palette(set, false)
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}
    p.Legend.Add(strings.Trim(asLabel[set], " pras"), l)
  }

  savePlot(p, "Anti-Stokes Background Subtracted")

  // Stokes
  title = "Stokes"
  xlabel = "Frequency (GHz)"
  ylabel = "Spectral Density (nV)"
  legend = "Pump"
  xrange = []float64{1, 1.36}
  yrange = []float64{0, 17.5}
  xtick = []float64{1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}
  ytick = []float64{0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5}
  xtickLabels = []string{"1", "", "1.1", "", "1.2", "", "1.3", "", "1.4"}
  ytickLabels = []string{"0", "", "5", "", "10", "", "15", ""}

  p = prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabels, ytickLabels,
  )

  for _, set := range sets {

    sPts := buildData(s[set])

    // Make a scatter plotter and set its style.
    plotSet, err := plotter.NewScatter(sPts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotSet.GlyphStyle.Color = palette(set, false)
    plotSet.GlyphStyle.Radius = vg.Points(3)
    plotSet.Shape = draw.CircleGlyph{}

    p.Add(plotSet)

    // Legend
    l, err := plotter.NewScatter(sPts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    l.GlyphStyle.Color = palette(set, false)
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}

    p.Legend.Add(strings.Trim(sLabel[set], " rs"), l)
  }

  savePlot(p, "Stokes Background Subtracted")
}

func goPlotasFits(
  sets []int,
  as, fits, widthLines [][][]float64,
  labels []string,
  widths, notes []float64,
  temp bool,
  sample string,
) {

  title := sample + "Anti-Stokes"
  xlabel := "Frequency (GHz)"
  ylabel := "Spectral Density (nV)"
  legend := "Power"

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("fits", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  for i, set := range sets {

    pts := buildData(as[set])
    fit := buildData(fits[i])
    wid := buildData(widthLines[i])

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotPts.GlyphStyle.Color = palette(set, false)
    plotPts.GlyphStyle.Radius = vg.Points(3)
    plotPts.Shape = draw.CircleGlyph{}

    // Plot fit
    plotFit, err := plotter.NewLine(fit)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotFit.LineStyle.Color = palette(set, true)
    plotFit.LineStyle.Width = vg.Points(3)

    // Width lines
    plotWid, err := plotter.NewLine(wid)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotWid.LineStyle.Color = palette(set, true)
    plotWid.LineStyle.Width = vg.Points(4)
    plotWid.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

    // Add set plots to p
    p.Add(plotPts, plotFit, plotWid)

    // Legend
    l, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    l.GlyphStyle.Color = palette(set, true)
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}
    power := strings.Trim(labels[set], " pras")
    if temp {
      temp := strconv.FormatFloat(notes[set], 'f', -1, 64)
      p.Legend.Add(power + " @" + temp + "K", l)
    } else {
      p.Legend.Add(power, l)
    }

  }

  savePlot(p, "Anti-Stokes w Fits")
}

func bin(
  as, s [][][]float64,
  binMHz float64,
) (
  [][][]float64, [][][]float64,
) {

  var asBinned, sBinned [][][]float64
  // copy as, s first [] to asBinned, sBinned w/append

  binGHz := binMHz/1000
  nBins := (as[0][0][len(as[0][0]) - 1] - as[0][0][0])/binGHz
  bins := []float64{}
  bound := as[0][0][0]

  for i := 0; i < int(nBins) + 1; i++ {
    bound += binGHz
    bins = append(bins, bound)

    var binFreqs []float64
    for _, f := range as[0][0] {
      if f < bound {
        binFreqs = append(binFreqs, f)
      }
    }
    // Want to average *signal* not frequencies
    // Frequency for each bin should be center of the bin
    // Create asSlice[0 = f, 1 = sig][vals], sSlice[0 = f, 1 = sig][vals]
    // append asBinned[set], sBinned[set] with asSlice, sSlice respectively  
    asBinned[0][0][i] = avg(binFreqs)
  }

  fmt.Printf("%.2f", asBinned[0][0])

  return asBinned, sBinned
}

func avg(
  toAvg []float64,
) (
  float64,
) {

  sum := 0.
  for _, v := range toAvg {
    sum += v
  }
  return sum/float64(len(toAvg))
}

func goPlotasPowerVsWid(
  sets []int,
  labels []string,
  notes, widths []float64,
  temp bool,
  sample string,
) {

  title := "Anti-Stokes Pump Power vs Widths of Fits"
  xlabel := "Pump Power (mW)"
  ylabel := "Full Width Half Max (MHz)"
  legend := ""

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("pow vs wid", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  for i, set := range sets {

    pts := make(plotter.XYs, 1)

    power := strings.Trim(labels[set], " mW pras")
    if pwr, err := strconv.ParseFloat(power, 64); err == nil {
      pts[0].X = pwr
    } else {
      fmt.Println(err)
      os.Exit(1)
    }
    pts[0].Y = widths[i]

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotPts.GlyphStyle.Color = palette(set, true)
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
      fmt.Println(err)
      os.Exit(1)
    }

    vDash.LineStyle.Color = palette(set, true)
    vDash.LineStyle.Width = vg.Points(4)
    vDash.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

    // Horizontal
    h[0].X = -15
    h[0].Y = pts[0].Y
    h[1].X = pts[0].X
    h[1].Y = pts[0].Y

    hDash, err := plotter.NewLine(h)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
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
      p.Legend.Add(power + " mW", plotPts)
    }
  }

  savePlot(p, "as Pow vs Wid")
}

func goPlotsFits(
  sets []int,
  s, fits, widthLines [][][]float64,
  labels []string,
  widths, notes []float64,
  temp bool,
  sample string,
) {

  title := "Stokes"
  xlabel := "Frequency (GHz)"
  ylabel := "Spectral Density (nV)"
  legend := "Power"

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("fits", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  for i, set := range sets {

    pts := buildData(s[set])
    fit := buildData(fits[i])
    wid := buildData(widthLines[i])

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotPts.GlyphStyle.Color = palette(set, false)
    plotPts.GlyphStyle.Radius = vg.Points(3)
    plotPts.Shape = draw.CircleGlyph{}

    // Plot fit
    plotFit, err := plotter.NewLine(fit)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotFit.LineStyle.Color = palette(set, true)
    plotFit.LineStyle.Width = vg.Points(3)

    // Width lines
    plotWid, err := plotter.NewLine(wid)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotWid.LineStyle.Color = palette(set, true)
    plotWid.LineStyle.Width = vg.Points(4)
    plotWid.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

    // Add set plots to p
    p.Add(plotPts, plotFit, plotWid)

    // Legend
    l, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    l.GlyphStyle.Color = palette(set, true)
    l.GlyphStyle.Radius = vg.Points(6)
    l.Shape = draw.CircleGlyph{}
    power := strings.Trim(labels[set], " prs")
    if temp {
      temperature := strconv.FormatFloat(notes[set], 'f', -1, 64)
      p.Legend.Add(power + " @" + temperature + "K", l)
    } else {
      p.Legend.Add(power, l)
    }
  }

  savePlot(p, "Stokes w Fits")
}

func goPlotsPowerVsWid(
  sets []int,
  labels []string,
  notes, widths []float64,
  temp bool,
  sample string,
) {

  title := "Stokes Pump Power vs Widths of Fits"
  xlabel := "Pump Power (mW)"
  ylabel := "Full Width Half Max (MHz)"
  legend := ""

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("pow vs wid", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  for i, set := range sets {

    pts := make(plotter.XYs, 1)

    power := strings.Trim(labels[set], " mW prs")
    if pwr, err := strconv.ParseFloat(power, 64); err == nil {
      pts[0].X = pwr
    } else {
      panic(err)
    }
    pts[0].Y = widths[i]

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotPts.GlyphStyle.Color = palette(set, true)
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
      fmt.Println(err)
      os.Exit(1)
    }

    vDash.LineStyle.Color = palette(set, true)
    vDash.LineStyle.Width = vg.Points(4)
    vDash.LineStyle.Dashes = []vg.Length{vg.Points(15), vg.Points(5)}

    // Horizontal
    h[0].X = -15
    h[0].Y = pts[0].Y
    h[1].X = pts[0].X
    h[1].Y = pts[0].Y

    hDash, err := plotter.NewLine(h)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
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

  savePlot(p, "s Pow vs Wid")
}

func goPlotHeightRatios(
  sets []int,
  heightRatios, powers []float64,
  labels []string,
  sample string,
) {

  title := "Height Ratios vs Power"
  xlabel := "Pump Power (mW)"
  ylabel := "Stokes/Anti-Stokes Heights"
  legend := ""

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("height ratios", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  // Linear fit line

  // Fit parameter guesses
  m := 5.
  b := .0125

  f := func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for i, set := range sets {

      x = powers[set]
      y := heightRatios[i]

      dst[i] = m * x + b - y
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

  results, err := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  m, b = results.X[0], results.X[1]

  var yFit []float64
  var xFit []float64

  // Create function according to solved fit parameters
  for _, set := range sets {
    var x float64

    x = powers[set]

    xFit = append(xFit, x)
    yFit = append(yFit, m * x + b)
  }

  fit := buildData([][]float64{xFit, yFit})

  // Plot fit
  plotFit, err := plotter.NewLine(fit)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p.Add(plotFit)

  plotFit.LineStyle.Color = color.RGBA{R: 127, G: 127, B: 127, A: 255}
  plotFit.LineStyle.Width = vg.Points(3)

  for i, set := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = powers[set]
    pts[0].Y = heightRatios[i]

    // Plot points
    plotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    plotPts.GlyphStyle.Color = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    plotPts.GlyphStyle.Radius = vg.Points(6)
    plotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(plotPts)
    //p.Legend.Add(strings.Trim(labels[set], " prs"), plotPts)
  }

  savePlot(p, "height ratios")
}

func goPlotLinewidths(
  sets []int,
  asLinewidths, sLinewidths, asPowers, sPowers []float64,
  labels []string,
  sample string,
) {

  title := "Linewidths vs Power"
  xlabel := "Pump Power (mW)"
  ylabel := "Full Width Half Max (MHz)"
  legend := ""

  xrange, yrange, xtick, ytick, xtickLabel, ytickLabel, err := axes("linewidths", sample)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  p := prepPlot(
    title, xlabel, ylabel, legend,
    xrange, yrange, xtick, ytick,
    xtickLabel, ytickLabel,
  )

  // as linear fit

  // Fit parameter guesses
  m := 1.
  b := 100.

  f := func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for i, set := range sets {

      x = asPowers[set]
      y := asLinewidths[i]

      dst[i] = m * x + b - y
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

  results, err := lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  m, b = results.X[0], results.X[1]

  var asyFit []float64
  var asxFit []float64

  // Create function according to solved fit parameters
  for _, set := range sets {
    var x float64

    x = asPowers[set]

    asxFit = append(asxFit, x)
    asyFit = append(asyFit, m * x + b)
  }

  asfit := buildData([][]float64{asxFit, asyFit})

  // Plot as fit
  asPlotFit, err := plotter.NewLine(asfit)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  asPlotFit.LineStyle.Color = color.RGBA{R: 99, G: 124, B: 198, A: 255}
  asPlotFit.LineStyle.Width = vg.Points(3)

  // s linear fit

  f = func(dst, guess []float64) {

    var x float64
    m, b := guess[0], guess[1]

    for i, set := range sets {

      x = sPowers[set]
      y := sLinewidths[i]

      dst[i] = m * x + b - y
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

  results, err = lm.LM(toBeSolved, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  m, b = results.X[0], results.X[1]

  var syFit []float64
  var sxFit []float64

  // Create function according to solved fit parameters
  for _, set := range sets {
    var x float64

    x = sPowers[set]

    sxFit = append(sxFit, x)
    syFit = append(syFit, m * x + b)
  }

  sfit := buildData([][]float64{sxFit, syFit})

  // Plot as fit
  sPlotFit, err := plotter.NewLine(sfit)
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  sPlotFit.LineStyle.Color = color.RGBA{R: 201, G: 104, B: 146, A: 255}
  sPlotFit.LineStyle.Width = vg.Points(3)

  p.Add(asPlotFit, sPlotFit)
  p.Legend.Add("Anti-Stokes", asPlotFit)
  p.Legend.Add("Stokes", sPlotFit)

  // as points
  for i, set := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = asPowers[set]
    pts[0].Y = asLinewidths[i]

    // Plot points
    asPlotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    asPlotPts.GlyphStyle.Color = color.RGBA{R: 99, G: 124, B: 198, A: 255}
    asPlotPts.GlyphStyle.Radius = vg.Points(6)
    asPlotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(asPlotPts)
  }

  // s points
  for i, set := range sets {

    pts := make(plotter.XYs, 1)

    pts[0].X = sPowers[set]
    pts[0].Y = sLinewidths[i]

    // Plot points
    sPlotPts, err := plotter.NewScatter(pts)
    if err != nil {
      fmt.Println(err)
      os.Exit(1)
    }

    sPlotPts.GlyphStyle.Color = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    sPlotPts.GlyphStyle.Radius = vg.Points(6)
    sPlotPts.Shape = draw.CircleGlyph{}

    // Add set plots to p
    p.Add(sPlotPts)
  }

  savePlot(p, "linewidths")
}

func prepPlot(
  title, xlabel, ylabel, legend string,
  xrange, yrange, xtick, ytick []float64,
  xtickLabels, ytickLabels []string,
) (
  *plot.Plot,
) {

  p := plot.New()
  p.BackgroundColor = color.RGBA{A:0}
  p.Title.Text = title
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 50
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = xlabel
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 36
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Min = xrange[0]
  p.X.Max = xrange[1]
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 36
  p.X.Tick.Label.Font.Variant = "Sans"

  xticks := []plot.Tick{}
  for i, v := range xtick {
    xticks = append(xticks, plot.Tick{Value: v, Label: xtickLabels[i]})
  }

  p.X.Tick.Marker = plot.ConstantTicks(xticks)
  p.X.Padding = vg.Points(-12.5)

  p.Y.Label.Text = ylabel
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 36
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Min = yrange[0]
  p.Y.Max = yrange[1]
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 36
  p.Y.Tick.Label.Font.Variant = "Sans"

  yticks := []plot.Tick{}
  for i, v := range ytick {
    yticks = append(yticks, plot.Tick{Value: v, Label: ytickLabels[i]})
  }

  p.Y.Tick.Marker = plot.ConstantTicks(yticks)
  p.Y.Padding = vg.Points(-4.75)

  p.Legend.TextStyle.Font.Size = 36
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)
  p.Legend.Add(legend)

  return p
}

func palette(
  brush int,
  dark bool,
) (
  color.RGBA,
) {

  if dark {
    darkColor := make([]color.RGBA, 16)
    darkColor[0] = color.RGBA{R: 27, G: 170, B: 139, A: 255}
    darkColor[4] = color.RGBA{R: 201, G: 104, B: 146, A: 255}
    darkColor[8] = color.RGBA{R: 99, G: 124, B: 198, A: 255}
    darkColor[12] = color.RGBA{R: 183, G: 139, B: 89, A: 255}
    darkColor[15] = color.RGBA{R: 18, G: 102, B: 99, A: 255}
    darkColor[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
    darkColor[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
    darkColor[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
    darkColor[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
    darkColor[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
    darkColor[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
    darkColor[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
    darkColor[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
    darkColor[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
    darkColor[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
    darkColor[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

    return darkColor[brush]
  }

  col := make([]color.RGBA, 16)
  col[0] = color.RGBA{R: 31, G: 211, B: 172, A: 255}
  col[4] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
  col[8] = color.RGBA{R: 122, G: 156, B: 255, A: 255}
  col[12] = color.RGBA{R: 255, G: 193, B: 122, A: 255}
  col[15] = color.RGBA{R: 27, G: 150, B: 146, A: 255}
  col[1] = color.RGBA{R: 188, G: 117, B: 255, A: 255}
  col[5] = color.RGBA{R: 234, G: 156, B: 172, A: 255}
  col[6] = color.RGBA{R: 1, G: 56, B: 84, A: 255}
  col[7] = color.RGBA{R: 46, G: 140, B: 60, A: 255}
  col[2] = color.RGBA{R: 140, G: 46, B: 49, A: 255}
  col[9] = color.RGBA{R: 122, G: 41, B: 104, A: 255}
  col[10] = color.RGBA{R: 41, G: 122, B: 100, A: 255}
  col[11] = color.RGBA{R: 122, G: 90, B: 41, A: 255}
  col[3] = color.RGBA{R: 91, G: 22, B: 22, A: 255}
  col[13] = color.RGBA{R: 22, G: 44, B: 91, A: 255}
  col[14] = color.RGBA{R: 59, G: 17, B: 66, A: 255}

  return col[brush]
}

func savePlot(
  p *plot.Plot, name string,
) {

  date := time.Now()

  // Make current date folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02"), 0755); err != nil {
      fmt.Println(err)
      os.Exit(1)
    }
  }

  // Make current time folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05"), 0755); err != nil {
      fmt.Println(err)
      os.Exit(1)
    }
  }

  path := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05") + "/" + name

  if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".png"); err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".svg"); err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, path + ".pdf"); err != nil {
    fmt.Println(err)
    os.Exit(1)
  }
}

func normalizeFit(
  fit []float64,
) (
  []float64,
) {

  var shift float64 = (fit[0] + fit[599])/2

  for i := range fit {
    fit[i] = fit[i] - shift
  }
  return fit
}

func writeLog(
  log []string,
) {

  date := time.Now()

  // Make current date folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02"), 0755); err != nil {
      fmt.Println(err)
      os.Exit(1)
    }
  }

  // Make current time folder if it doesn't already exist
  if _, err := os.Stat("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")); os.IsNotExist(err) {
    if err := os.Mkdir("plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05"), 0755); err != nil {
      fmt.Println(err)
      os.Exit(1)
    }
  }

  path := "plots/" + date.Format("2006-Jan-02") + "/" + date.Format("15:04:05")

  txt, err := os.Create(path + "/log.txt")
  if err != nil {
    fmt.Println(err)
    os.Exit(1)
  }

  w := bufio.NewWriter(txt)
  defer w.Flush()
  for _, line := range log {
    if _, err := w.WriteString(line); err != nil {
      fmt.Println(err)
      os.Exit(1)
    }
  }
}
