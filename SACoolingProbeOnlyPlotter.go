package main

import (
  "github.com/Arafatk/glot"
  "encoding/csv"
  "bufio"
  "fmt"
  "os"
  "io"
  "strconv"
  "strings"
  "math"
  "gonum.org/v1/gonum/floats"
	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/gonum/diff/fd"
	"gonum.org/v1/gonum/optimize"
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

  setsToPlotSubtracted := []int{3}
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

  setsToPlotSubtractedGrouped := []int{}
  plotSubtractedGrouped(
    setsToPlotSubtractedGrouped,
    s, rsLabel,
    as, rasLabel,
  )
/*
  // Lorentz Fit
  asamp := []float64{0, 0.01, 0.03, 0.1, 0.5, 0.7, 0.1}
  amp := 0.1
  width := 0.1
  center := 2.27

  //fitParams := []float64{amp, width, center}

  as1LorentzJac := NumJac{Func: as1Lorentz}
  s1LorentzJac := NumJac{Func: s1Lorentz}
  as2LorentzJac := NumJac{Func: as2Lorentz}
  s2LorentzJac := NumJac{Func: s2Lorentz}
  as3LorentzJac := NumJac{Func: as3Lorentz}
  s3LorentzJac := NumJac{Func: s3Lorentz}
  as4LorentzJac := NumJac{Func: as4Lorentz}
  s4LorentzJac := NumJac{Func: s4Lorentz}
  as5LorentzJac := NumJac{Func: as5Lorentz}
  s5LorentzJac := NumJac{Func: s5Lorentz}
  as6LorentzJac := NumJac{Func: as6Lorentz}
  s6LorentzJac := NumJac{Func: s6Lorentz}

  as1LorentzProb := LMProblem{
	  Dim:        3,
 	  Size:       len(as1data[0]),
 	  Func:       as1Lorentz,
 	  Jac:        as1LorentzJac.Jac,
 	  InitParams: []float64{asamp[1], width, center},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }
  s1LorentzProb := LMProblem{
	  Dim:        3,
 	  Size:       len(s1data[0]),
 	  Func:       s1Lorentz,
 	  Jac:        s1LorentzJac.Jac,
 	  InitParams: []float64{amp, width, center},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }
  as2LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(as2data[0]),
    Func:       as2Lorentz,
    Jac:        as2LorentzJac.Jac,
    InitParams: []float64{asamp[2], width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  s2LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(s2data[0]),
    Func:       s2Lorentz,
    Jac:        s2LorentzJac.Jac,
    InitParams: []float64{amp, width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  as3LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(as3data[0]),
    Func:       as3Lorentz,
    Jac:        as3LorentzJac.Jac,
    InitParams: []float64{asamp[3], width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  s3LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(s3data[0]),
    Func:       s3Lorentz,
    Jac:        s3LorentzJac.Jac,
    InitParams: []float64{amp, width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  as4LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(as4data[0]),
    Func:       as4Lorentz,
    Jac:        as4LorentzJac.Jac,
    InitParams: []float64{asamp[4], width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  s4LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(s4data[0]),
    Func:       s4Lorentz,
    Jac:        s4LorentzJac.Jac,
    InitParams: []float64{amp, width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  as5LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(as5data[0]),
    Func:       as5Lorentz,
    Jac:        as5LorentzJac.Jac,
    InitParams: []float64{asamp[5], width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  s5LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(s5data[0]),
    Func:       s5Lorentz,
    Jac:        s5LorentzJac.Jac,
    InitParams: []float64{amp, width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  as6LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(as6data[0]),
    Func:       as6Lorentz,
    Jac:        as6LorentzJac.Jac,
    InitParams: []float64{asamp[6], width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }
  s6LorentzProb := LMProblem{
    Dim:        3,
    Size:       len(s6data[0]),
    Func:       s6Lorentz,
    Jac:        s6LorentzJac.Jac,
    InitParams: []float64{amp, width, center},
    Tau:        1e-6,
    Eps1:       1e-8,
    Eps2:       1e-8,
  }

  as1LorentzResults, _ := LM(as1LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s1LorentzResults, _ := LM(s1LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  as2LorentzResults, _ := LM(as2LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s2LorentzResults, _ := LM(s2LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  as3LorentzResults, _ := LM(as3LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s3LorentzResults, _ := LM(s3LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  as4LorentzResults, _ := LM(as4LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s4LorentzResults, _ := LM(s4LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  as5LorentzResults, _ := LM(as5LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s5LorentzResults, _ := LM(s5LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  as6LorentzResults, _ := LM(as6LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s6LorentzResults, _ := LM(s6LorentzProb, &Settings{Iterations: 100, ObjectiveTol: 1e-16})

  fmt.Println(as1LorentzResults.X)
  fmt.Println(s1LorentzResults.X)
  fmt.Println(as2LorentzResults.X)
  fmt.Println(s2LorentzResults.X)
  fmt.Println(as3LorentzResults.X)
  fmt.Println(s3LorentzResults.X)
  fmt.Println(as4LorentzResults.X)
  fmt.Println(s4LorentzResults.X)
  fmt.Println(as5LorentzResults.X)
  fmt.Println(s5LorentzResults.X)
  fmt.Println(as6LorentzResults.X)
  fmt.Println(s6LorentzResults.X)

  var as1yfit []float64
  var s1yfit []float64
  var as2yfit []float64
  var s2yfit []float64
  var as3yfit []float64
  var s3yfit []float64
  var as4yfit []float64
  var s4yfit []float64
  var as5yfit []float64
  var s5yfit []float64
  var as6yfit []float64
  var s6yfit []float64

  for i := 0; i < 600; i++ {
    // (amp*wid^2/((x-cen)^2+wid^2))
    as1yfit = append(as1yfit, as1LorentzResults.X[0] * math.Pow(as1LorentzResults.X[1], 2) / (math.Pow(as1data[0][i] - as1LorentzResults.X[2], 2) + math.Pow(as1LorentzResults.X[1], 2)))
    s1yfit = append(s1yfit, s1LorentzResults.X[0] * math.Pow(s1LorentzResults.X[1], 2) / (math.Pow(s1data[0][i] - s1LorentzResults.X[2], 2) + math.Pow(s1LorentzResults.X[1], 2)))
    as2yfit = append(as2yfit, as2LorentzResults.X[0] * math.Pow(as2LorentzResults.X[1], 2) / (math.Pow(as2data[0][i] - as2LorentzResults.X[2], 2) + math.Pow(as2LorentzResults.X[1], 2)))
    s2yfit = append(s2yfit, s2LorentzResults.X[0] * math.Pow(s2LorentzResults.X[1], 2) / (math.Pow(s2data[0][i] - s2LorentzResults.X[2], 2) + math.Pow(s2LorentzResults.X[1], 2)))
    as3yfit = append(as3yfit, as3LorentzResults.X[0] * math.Pow(as3LorentzResults.X[1], 2) / (math.Pow(as3data[0][i] - as3LorentzResults.X[2], 2) + math.Pow(as3LorentzResults.X[1], 2)))
    s3yfit = append(s3yfit, s3LorentzResults.X[0] * math.Pow(s3LorentzResults.X[1], 2) / (math.Pow(s3data[0][i] - s3LorentzResults.X[2], 2) + math.Pow(s3LorentzResults.X[1], 2)))
    as4yfit = append(as4yfit, as4LorentzResults.X[0] * math.Pow(as4LorentzResults.X[1], 2) / (math.Pow(as4data[0][i] - as4LorentzResults.X[2], 2) + math.Pow(as4LorentzResults.X[1], 2)))
    s4yfit = append(s4yfit, s4LorentzResults.X[0] * math.Pow(s4LorentzResults.X[1], 2) / (math.Pow(s4data[0][i] - s4LorentzResults.X[2], 2) + math.Pow(s4LorentzResults.X[1], 2)))
    as5yfit = append(as5yfit, as5LorentzResults.X[0] * math.Pow(as5LorentzResults.X[1], 2) / (math.Pow(as5data[0][i] - as5LorentzResults.X[2], 2) + math.Pow(as5LorentzResults.X[1], 2)))
    s5yfit = append(s5yfit, s5LorentzResults.X[0] * math.Pow(s5LorentzResults.X[1], 2) / (math.Pow(s5data[0][i] - s5LorentzResults.X[2], 2) + math.Pow(s5LorentzResults.X[1], 2)))
    as6yfit = append(as6yfit, as6LorentzResults.X[0] * math.Pow(as6LorentzResults.X[1], 2) / (math.Pow(as6data[0][i] - as6LorentzResults.X[2], 2) + math.Pow(as6LorentzResults.X[1], 2)))
    s6yfit = append(s6yfit, s6LorentzResults.X[0] * math.Pow(s6LorentzResults.X[1], 2) / (math.Pow(s6data[0][i] - s6LorentzResults.X[2], 2) + math.Pow(s6LorentzResults.X[1], 2)))
  }

  normalizeFit(as1yfit)
  normalizeFit(s1yfit)
  normalizeFit(as2yfit)
  normalizeFit(s2yfit)
  normalizeFit(as3yfit)
  normalizeFit(s3yfit)
  normalizeFit(as4yfit)
  normalizeFit(s4yfit)
  normalizeFit(as5yfit)
  normalizeFit(s5yfit)
  normalizeFit(as6yfit)
  normalizeFit(s6yfit)

  as1Fit := [][]float64{as1data[0], as1yfit}
  s1Fit := [][]float64{s1data[0], s1yfit}
  as2Fit := [][]float64{as2data[0], as2yfit}
  s2Fit := [][]float64{s2data[0], s2yfit}
  as3Fit := [][]float64{as3data[0], as3yfit}
  s3Fit := [][]float64{s3data[0], s3yfit}
  as4Fit := [][]float64{as4data[0], as4yfit}
  s4Fit := [][]float64{s4data[0], s4yfit}
  as5Fit := [][]float64{as5data[0], as5yfit}
  s5Fit := [][]float64{s5data[0], s5yfit}
  as6Fit := [][]float64{as6data[0], as6yfit}
  s6Fit := [][]float64{s6data[0], s6yfit}

  plotDataWithFit(
    as1, label[2],
    as1Fit, label[2],
    s1, label[4],
    s1Fit, label[4],
    as2, label[6],
    as2Fit, label[6],
    s2, label[8],
    s2Fit, label[8],
    as3, label[10],
    as3Fit, label[10],
    s3, label[12],
    s3Fit, label[12],
    as4, label[14],
    as4Fit, label[14],
    s4, label[16],
    s4Fit, label[16],
    as5, label[18],
    as5Fit, label[18],
    s5, label[20],
    s5Fit, label[20],
    as6, label[22],
    as6Fit, label[22],
    s6, label[24],
    s6Fit, label[24],
  )

  plotSeparateFits(
    as1Fit, label[2],
    as2Fit, label[6],
    as3Fit, label[10],
    as4Fit, label[14],
    as5Fit, label[18],
    as6Fit, label[22],
  )

  plotSeparateFits(
    s1Fit, label[4],
    s2Fit, label[8],
    s3Fit, label[12],
    s4Fit, label[16],
    s5Fit, label[20],
    s6Fit, label[24],
  )

  plotFits(
    as1Fit, label[2],
    s1Fit, label[4],
    as2Fit, label[6],
    s2Fit, label[8],
    as3Fit, label[10],
    s3Fit, label[12],
    as4Fit, label[14],
    s4Fit, label[16],
    as5Fit, label[18],
    s5Fit, label[20],
    as6Fit, label[22],
    s6Fit, label[24],
  )
  */
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

  // Background
  basLabel = append(basLabel, label[1])
  bsLabel = append(bsLabel, label[2])

  for i := 3; i < len(label); i += 2 {
    rasLabel = append(rasLabel, label[i])
  }
  for i := 4; i < len(label); i += 2 {
    rsLabel = append(rsLabel, label[i])
  }

  return rasLabel, basLabel, rsLabel, bsLabel
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

func plotSubtractedGrouped(
  sets []int,
  s [][][]float64, sLabel []string,
  as [][][]float64, asLabel []string,
  ) {

  // s
  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  for i := 0; i < len(sets); i++ {
    plot.AddPointGroup(strings.Trim(sLabel[sets[i]], " prs") + " s", "points", s[sets[i]])
  }

  // as
  dimensions = 2
  persist = true
  debug = false
  plot, _ = glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  for i := 0; i < len(sets); i++ {
    plot.AddPointGroup(strings.Trim(asLabel[sets[i]], " pras") + " as", "points", as[sets[i]])
  }
}

func normalizeFit(fit []float64) ([]float64) {

  var shift float64 = (fit[0] + fit[599])/2

  for i := 0; i < 600; i++ {
    fit[i] = fit[i] - shift
  }
  return fit
}

func plotDataWithFit(
  set1 [][]float64, label1 string,
  set2 [][]float64, label2 string,
  set3 [][]float64, label3 string,
  set4 [][]float64, label4 string,
  set5 [][]float64, label5 string,
  set6 [][]float64, label6 string,
  set7 [][]float64, label7 string,
  set8 [][]float64, label8 string,
  set9 [][]float64, label9 string,
  set10 [][]float64, label10 string,
  set11 [][]float64, label11 string,
  set12 [][]float64, label12 string,
  set13 [][]float64, label13 string,
  set14 [][]float64, label14 string,
  set15 [][]float64, label15 string,
  set16 [][]float64, label16 string,
  set17 [][]float64, label17 string,
  set18 [][]float64, label18 string,
  set19 [][]float64, label19 string,
  set20 [][]float64, label20 string,
  set21 [][]float64, label21 string,
  set22 [][]float64, label22 string,
  set23 [][]float64, label23 string,
  set24 [][]float64, label24 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted with Lorentzian Fit")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  plot.AddPointGroup(label1, "points", set1)
  plot.AddPointGroup(label2, "lines", set2)
  plot.AddPointGroup(label3, "points", set3)
  plot.AddPointGroup(label4, "lines", set4)
  plot.AddPointGroup(label5, "points", set5)
  plot.AddPointGroup(label6, "lines", set6)
  plot.AddPointGroup(label7, "points", set7)
  plot.AddPointGroup(label8, "lines", set8)
  plot.AddPointGroup(label9, "points", set9)
  plot.AddPointGroup(label10, "lines", set10)
  plot.AddPointGroup(label11, "points", set11)
  plot.AddPointGroup(label12, "lines", set12)
  plot.AddPointGroup(label13, "points", set13)
  plot.AddPointGroup(label14, "lines", set14)
  plot.AddPointGroup(label15, "points", set15)
  plot.AddPointGroup(label16, "lines", set16)
  plot.AddPointGroup(label17, "points", set17)
  plot.AddPointGroup(label18, "lines", set18)
  plot.AddPointGroup(label19, "points", set19)
  plot.AddPointGroup(label20, "lines", set20)
  plot.AddPointGroup(label21, "points", set21)
  plot.AddPointGroup(label22, "lines", set22)
  plot.AddPointGroup(label23, "points", set23)
  plot.AddPointGroup(label24, "lines", set24)
}

func plotSeparateFits(
set1 [][]float64, label1 string,
set2 [][]float64, label2 string,
set3 [][]float64, label3 string,
set4 [][]float64, label4 string,
set5 [][]float64, label5 string,
set6 [][]float64, label6 string,
) {

dimensions := 2
persist := true
debug := false
plot, _ := glot.NewPlot(dimensions, persist, debug)

plot.SetTitle("Separated")
plot.SetXLabel("Frequency (GHz)")
plot.SetYLabel("Signal (uV)")

plot.AddPointGroup(label1, "lines", set1)
plot.AddPointGroup(label2, "lines", set2)
plot.AddPointGroup(label3, "lines", set3)
plot.AddPointGroup(label4, "lines", set4)
plot.AddPointGroup(label5, "lines", set5)
plot.AddPointGroup(label6, "lines", set6)
}

func plotFits(
set1 [][]float64, label1 string,
set2 [][]float64, label2 string,
set3 [][]float64, label3 string,
set4 [][]float64, label4 string,
set5 [][]float64, label5 string,
set6 [][]float64, label6 string,
set7 [][]float64, label7 string,
set8 [][]float64, label8 string,
set9 [][]float64, label9 string,
set10 [][]float64, label10 string,
set11 [][]float64, label11 string,
set12 [][]float64, label12 string,
) {

dimensions := 2
persist := true
debug := false
plot, _ := glot.NewPlot(dimensions, persist, debug)

plot.SetTitle("Lorentz Fit")
plot.SetXLabel("Frequency (GHz)")
plot.SetYLabel("Signal (uV)")

plot.AddPointGroup(label1, "lines", set1)
plot.AddPointGroup(label2, "lines", set2)
plot.AddPointGroup(label3, "lines", set3)
plot.AddPointGroup(label4, "lines", set4)
plot.AddPointGroup(label5, "lines", set5)
plot.AddPointGroup(label6, "lines", set6)
plot.AddPointGroup(label7, "lines", set7)
plot.AddPointGroup(label8, "lines", set8)
plot.AddPointGroup(label9, "lines", set9)
plot.AddPointGroup(label10, "lines", set10)
plot.AddPointGroup(label11, "lines", set11)
plot.AddPointGroup(label12, "lines", set12)
}

/*
func lorentzFit(data [][]float64, fitParams []float64) ([][]float64) {

  lorentzJac := lm.NumJac{Func: set}

  lorentzProb := lm.LMProblem{
	  Dim:        3,
 	  Size:       len(data[0]),
 	  Func:       set,
 	  Jac:        lorentzJac.Jac,
 	  InitParams: []float64{fitParams[0], fitParams[1], fitParams[2]},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }

  lorentzResults, _ := lm.LM(lorentzProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  fmt.Println(lorentzResults.X)

  var yfit []float64

  for i := 0; i < 600; i++ {
    // (amp*wid^2/((x-cen)^2+wid^2))
    yfit = append(yfit, lorentzResults.X[0] * math.Pow(lorentzResults.X[1], 2) / (math.Pow(data[0][i] - lorentzResults.X[2], 2) + math.Pow(lorentzResults.X[1], 2)))
  }

  return [][]float64{data[0], yfit}
}
*/



/*
func as1Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[1])
  as1 := getData(file[2])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s1Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[3])
  s1 := getData(file[4])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func as2Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[5])
  as1 := getData(file[6])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s2Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[7])
  s1 := getData(file[8])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func as3Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[9])
  as1 := getData(file[10])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s3Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[11])
  s1 := getData(file[12])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func as4Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[13])
  as1 := getData(file[14])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s4Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[15])
  s1 := getData(file[16])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func as5Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[17])
  as1 := getData(file[18])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s5Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[19])
  s1 := getData(file[20])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func as6Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bas1 := getData(file[21])
  as1 := getData(file[22])
  as1Sub := subtractBackground(bas1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
func s6Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  bs1 := getData(file[23])
  s1 := getData(file[24])
  s1Sub := subtractBackground(bs1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
*/


/*
Package lm implements optimization routines for non-linear least squares problems
using the Levenberg-Marquardt method.

Given function f:Rn -> Rm, where m is the number of non-linear functions and n parameters,
the Levenberg-Marquardt method is used to seek a point X that minimizes F(x) = 0.5 * f.T * f.

The user supplies a non-linear function. The jacobian may also be supplied by the user or
approximated by finite differences.
*/

type Settings struct {
	// Iterations represents the maximum number of iterations. Defaults to 100.
	Iterations int

	// ObjectiveTol represents the value for the obejective function after which
	// the algorithm can stop. Defaults to 1e-16.
	ObjectiveTol float64
}

func defaultSettings(set *Settings) {
	set.Iterations = 100
	set.ObjectiveTol = 1e-16
}

type Result struct {
	X      []float64
	Status optimize.Status
}

// NumJac is used if the user doesn't wish to provide a fucnction that evaluates
// the jacobian matrix. NumJac provides a method Jac that computes the jacobian matrix
// by finite differences.
type NumJac struct {
	Func func(dst, param []float64)
}

func (nj *NumJac) Jac(dst *mat.Dense, param []float64) {
	fd.Jacobian(dst, nj.Func, param, &fd.JacobianSettings{
		Formula:    fd.Central,
		Concurrent: true,
	})
}

func maxDiagElem(m *mat.Dense) float64 {
	r, c := m.Dims()
	if r != c {
		panic("lm: matrix is not square")
	}
	maxElem := m.At(0, 0)
	for i := 1; i < r; i++ {
		if m.At(i, i) > maxElem {
			maxElem = m.At(i, i)
		}
	}
	return maxElem
}

func addToDiag(m *mat.Dense, v float64) {
	r, c := m.Dims()
	if r != c {
		panic("lm: matrix is not square")
	}
	for i := 0; i < r; i++ {
		m.Set(i, i, m.At(i, i)+v)
	}
}

func updateParams(dst []float64, params []float64, h *mat.VecDense) {
	if len(params) != h.Len() {
		panic("lm: lenghts don't match")
	}
	for i := 0; i < len(params); i++ {
		dst[i] = params[i] - h.At(i, 0)
	}
}

func calcRho(fParams []float64, fParamsNew []float64, h *mat.VecDense, grad *mat.VecDense, mu float64) float64 {
	rho := floats.Dot(fParams, fParams) - floats.Dot(fParamsNew, fParamsNew)
	tmpVec := mat.NewVecDense(h.Len(), nil)
	tmpVec.AddScaledVec(grad, mu, h)
	lDiff := mat.Dot(h, tmpVec)
	rho /= lDiff
	return rho
}

// LM is a function that solves non-linear least squares problems using the Levenberg-Marquardt
// Method.
//
// References:
//  - Madsen, Kaj, Hans Bruun Nielsen, and Ole Tingleff. "Methods for non-linear least squares
//    problems.", 2nd edition, 2004.
//  - Lourakis, Manolis. "A Brief Description of the Levenberg-Marquardt Algorithm Implemened
//    by levmar", 2005.
func LM(problem LMProblem, settings *Settings) (*Result, error) {
	var set Settings
	if settings != nil {
		set = *settings
	} else {
		defaultSettings(&set)
	}
	dim := problem.Dim
	if problem.Dim == 0 {
		panic("lm: problem dimension is 0")
	}
	size := problem.Size
	if problem.Size == 0 {
		panic("lm: problem size is 0")
	}
	status := optimize.NotTerminated

	dstFunc := make([]float64, size)
	dstFuncNew := make([]float64, size)
	dstJac := mat.NewDense(size, dim, nil)
	dstA := mat.NewDense(dim, dim, nil)
	dstGrad := mat.NewVecDense(dim, nil)
	dstH := mat.NewVecDense(dim, nil)
	nu := 2.0
	var mu float64
	found := false

	// The inital guess is the zero vector by default.
	parameters := make([]float64, dim)
	parametersNew := make([]float64, dim)
	if problem.InitParams != nil {
		copy(parameters, problem.InitParams)
	}

	// Initial evaluation of A = J.T * J and g = J.T * f.
	problem.Func(dstFunc, parameters)
	problem.Jac(dstJac, parameters)
	dstA.Mul(dstJac.T(), dstJac)
	dstGrad.MulVec(dstJac.T(), mat.NewVecDense(size, dstFunc))

	found = (mat.Norm(dstGrad, math.Inf(1)) <= problem.Eps1)
	mu = problem.Tau * maxDiagElem(dstA)

	for iter := 0; ; iter++ {
		if iter == set.Iterations {
			status = optimize.IterationLimit
			break
		}
		if found {
			status = optimize.StepConvergence
			break
		}

		// Solve (A + mu * I) * h_lm = g.
		addToDiag(dstA, mu)
		err := dstH.SolveVec(dstA, dstGrad)
		if err != nil {
			panic("singular")
		}

		// Return A to its original state for the next steps. This is done in order not to copy A.
		addToDiag(dstA, -mu)

		if mat.Norm(dstH, 2) <= (floats.Norm(parameters, 2)+problem.Eps2)*problem.Eps2 {
			found = true
		} else {
			updateParams(parametersNew, parameters, dstH)

			// Calculate rho = (F(x) - F(x_new)) / (L(0) - L(h_lm)), where
			// F = 0.5 * f.T * f, L = 0.5 * h_lm.T * (mu * h_lm - g).
			problem.Func(dstFuncNew, parametersNew)
			rho := calcRho(dstFunc, dstFuncNew, dstH, dstGrad, mu)

			if rho > 0 { // step is acceptable
				copy(parameters, parametersNew)
				problem.Func(dstFunc, parameters)
				problem.Jac(dstJac, parameters)
				dstA.Mul(dstJac.T(), dstJac)
				dstGrad.MulVec(dstJac.T(), mat.NewVecDense(size, dstFunc))
				found = (mat.Norm(dstGrad, math.Inf(1)) <= problem.Eps1) ||
					(0.5*floats.Dot(dstFunc, dstFunc) <= set.ObjectiveTol)
				mu = mu * math.Max(1.0/3.0, 1-math.Pow(2*rho-1, 3))
				nu = 2.0
			} else {
				mu *= nu
				nu *= 2.0
			}
		}
	}
	return &Result{
		X:      parameters,
		Status: status,
	}, nil
}

// LMProblem is used for running LM optimization. The objective function is
// F = 0.5 * f.T * f, where f:Rn -> Rm and m >= n.
type LMProblem struct {
	// Dim is the dimension of the parameters of the problem (n).
	Dim int
	// Size specifies the number of nonlinear functions (m).
	Size int
	// Func computes the function value at params.
	Func func(dst, param []float64)
	// Jac computes the jacobian matrix of Func.
	Jac func(dst *mat.Dense, param []float64)
	// InitParams stores the users inital guess. Defaults to the zero vector when nil.
	InitParams []float64
	// Tau scales the initial damping parameter.
	Tau float64
	// Eps1 is a stopping criterion for the gradient of F.
	Eps1 float64
	// Eps2 is a stopping criterion for the step size.
	Eps2 float64
}
