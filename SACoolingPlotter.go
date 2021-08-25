package main

import (
  "github.com/Arafatk/glot"
  "github.com/maorshutman/lm"
  "encoding/csv"
  "fmt"
  "os"
  "strconv"
  "strings"
  "math"
)

func main() {

  label, file := readMeta()

  background1 := getData(file[1])
  as1 := getData(file[2])
  //s1 := getData(file[3])
  //background2 := getData(file[4])
  //as2 := getData(file[5])
  //s2 := getData(file[6])

  plotRaw(
    background1, label[1],
    as1, label[2],
    //s1, label[3],
    //background2, label[4],
    //as2, label[5],
    //s2, label[6],
  )

  // Lorentz Fit
  amp := 0.4
  width := 0.1
  center := 2.25

  fitParams := []float64{amp, width, center}

  as1Fit := lorentzFit(subtractBackground(background1, as1), fitParams)
  //s1Fit := lorentzFit(subtractBackground(background1, s1), fitParams)

  plotData(
    as1, label[2],
    as1Fit, "Anti-Stokes Lorentz",
    //s1, label[3],
    //s1Fit, "Stokes Lorentz",
  )
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
    data = append(data, value[4])
  }

  return label, data
}

func getData(csvName string) ([][]float64) {

  // Read
  file, err := os.Open(csvName)
  if err != nil {
    fmt.Println(err)
  }

  reader := csv.NewReader(file)
  dataStr, err := reader.ReadAll()
  if err != nil {
    fmt.Println(err)
  }

  // Separate, Strip, & Transpose
  var frequencyStrT, signalStrT []string

  for i := 1; i < 601; i++ {
    frequencyStrT = append(frequencyStrT, strings.ReplaceAll(dataStr[i][0]," ",""))
    signalStrT = append(signalStrT, strings.ReplaceAll(dataStr[i][2]," ",""))
  }

  // Convert
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

  return [][]float64{frequency, signal}
}

func plotRaw(
  set1 [][]float64, label1 string,
  //set2 [][]float64, label2 string,
  set3 [][]float64, label3 string,
  //set4 [][]float64, label4 string,
  //set5 [][]float64, label5 string,
  //set6 [][]float64, label6 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Probe Only - Full Power - Big RF Amp w/Spectrum Analyzer")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (dBm)")

  plot.AddPointGroup(label1, "points", set1)
  //plot.AddPointGroup(label2, "points", set2)
  plot.AddPointGroup(label3, "points", set3)
}

func subtractBackground(b [][]float64, s [][]float64) ([][]float64) {

  for i := 0; i < 600; i++ {
    s[1][i] = s[1][i] - b[1][i]
  }

  return s
}

func plotData(
  set1 [][]float64, label1 string,
  //set2 [][]float64, label2 string,
  set3 [][]float64, label3 string,
  //set4 [][]float64, label4 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Probe Only - Full Power - Big RF Amp w/Spectrum Analyzer")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (dBm)")

  plot.AddPointGroup(label1, "points", set1)
  //plot.AddPointGroup(label2, "points", set2)
  plot.AddPointGroup(label3, "lines", set3)
  //plot.AddPointGroup(label4, "points", set4)
}

func lorentzFit(data [][]float64, fitParams []float64) ([][]float64) {

  lorentzJac := lm.NumJac{Func: lorentz}

  lorentzProb := lm.LMProblem{
	  Dim:        3,
 	  Size:       len(data[0]),
 	  Func:       lorentz,
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

  //plot, _ := glot.NewPlot(2, true, false)
  //plot.AddPointGroup("Lorentz Fit", "lines", fit)
}

func lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  background1 := getData(file[1])
  as1 := getData(file[2])
  data := subtractBackground(background1, as1)

  for i := 0; i < 600; i++ {
    x := data[0][i]
    y := data[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
