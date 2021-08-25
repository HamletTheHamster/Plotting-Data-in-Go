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
  s1 := getData(file[3])
  //background2 := getData(file[4])
  //as2 := getData(file[5])
  //s2 := getData(file[6])

  plotRaw(
    background1, label[1],
    as1, label[2],
    s1, label[3],
    //background2, label[4],
    //as2, label[5],
    //s2, label[6],
  )

  as1data := subtractBackground(background1, as1)
  s1data := subtractBackground(background1, s1)

  plotSubtracted(
    as1data, label[2],
    s1data, label[3],
  )

  // Lorentz Fit
  amp := 0.4
  width := 0.1
  center := 2.25

  //fitParams := []float64{amp, width, center}

  as1LorentzJac := lm.NumJac{Func: as1Lorentz}
  s1LorentzJac := lm.NumJac{Func: s1Lorentz}

  as1LorentzProb := lm.LMProblem{
	  Dim:        3,
 	  Size:       len(as1data[0]),
 	  Func:       as1Lorentz,
 	  Jac:        as1LorentzJac.Jac,
 	  InitParams: []float64{amp, width, center},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }
  s1LorentzProb := lm.LMProblem{
	  Dim:        3,
 	  Size:       len(s1data[0]),
 	  Func:       s1Lorentz,
 	  Jac:        s1LorentzJac.Jac,
 	  InitParams: []float64{amp, width, center},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }

  as1LorentzResults, _ := lm.LM(as1LorentzProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})
  s1LorentzResults, _ := lm.LM(s1LorentzProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  fmt.Println(as1LorentzResults.X)
  fmt.Println(s1LorentzResults.X)

  var as1yfit []float64
  var s1yfit []float64

  for i := 0; i < 600; i++ {
    // (amp*wid^2/((x-cen)^2+wid^2))
    as1yfit = append(as1yfit, as1LorentzResults.X[0] * math.Pow(as1LorentzResults.X[1], 2) / (math.Pow(as1data[0][i] - as1LorentzResults.X[2], 2) + math.Pow(as1LorentzResults.X[1], 2)))
    s1yfit = append(s1yfit, s1LorentzResults.X[0] * math.Pow(s1LorentzResults.X[1], 2) / (math.Pow(s1data[0][i] - s1LorentzResults.X[2], 2) + math.Pow(s1LorentzResults.X[1], 2)))
  }

  as1Fit := [][]float64{as1data[0], as1yfit}
  s1Fit := [][]float64{s1data[0], s1yfit}

  //as1Fit := lorentzFit(subtractBackground(background1, as1), fitParams)
  //s1Fit := lorentzFit(subtractBackground(background1, s1), fitParams)

  plotDataWithFit(
    as1, label[2],
    as1Fit, "Anti-Stokes Lorentz",
    s1, label[3],
    s1Fit, "Stokes Lorentz",
  )

  plotFits(
    as1Fit, "Anti-Stokes",
    s1Fit, "Stokes",
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
  set2 [][]float64, label2 string,
  set3 [][]float64, label3 string,
  //set4 [][]float64, label4 string,
  //set5 [][]float64, label5 string,
  //set6 [][]float64, label6 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Probe Only - Raw")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (dBm)")

  plot.AddPointGroup(label1, "points", set1)
  plot.AddPointGroup(label2, "points", set2)
  plot.AddPointGroup(label3, "points", set3)
}

func subtractBackground(b [][]float64, s [][]float64) ([][]float64) {

  for i := 0; i < 600; i++ {
    s[1][i] = s[1][i] - b[1][i]
  }

  return s
}

func plotSubtracted(
  set1 [][]float64, label1 string,
  set2 [][]float64, label2 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Probe Only - Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (dBm)")

  plot.AddPointGroup(label1, "points", set1)
  plot.AddPointGroup(label2, "points", set2)
}

func plotDataWithFit(
  set1 [][]float64, label1 string,
  set2 [][]float64, label2 string,
  set3 [][]float64, label3 string,
  set4 [][]float64, label4 string,
  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Probe Only - Background Subtracted with Lorentzian Fit")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (dBm)")

  plot.AddPointGroup(label1, "points", set1)
  plot.AddPointGroup(label2, "lines", set2)
  plot.AddPointGroup(label3, "points", set3)
  plot.AddPointGroup(label4, "lines", set4)
}

func plotFits(
set1 [][]float64, label1 string,
set2 [][]float64, label2 string,
) {

dimensions := 2
persist := true
debug := false
plot, _ := glot.NewPlot(dimensions, persist, debug)

plot.SetTitle("Probe Only - Lorentz Fit")
plot.SetXLabel("Frequency (GHz)")
plot.SetYLabel("Signal (dBm)")

plot.AddPointGroup(label1, "lines", set1)
plot.AddPointGroup(label2, "lines", set2)
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
func as1Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  background1 := getData(file[1])
  as1 := getData(file[2])
  as1Sub := subtractBackground(background1, as1)

  for i := 0; i < 600; i++ {
    x := as1Sub[0][i]
    y := as1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}

func s1Lorentz(dst, fitParams []float64) {

  _, file := readMeta()
  background1 := getData(file[1])
  s1 := getData(file[3])
  s1Sub := subtractBackground(background1, s1)

  for i := 0; i < 600; i++ {
    x := s1Sub[0][i]
    y := s1Sub[1][i]
    //fitParams = {amp, width, center}
    dst[i] = fitParams[0] * math.Pow(fitParams[1], 2) / (math.Pow(x - fitParams[2], 2) + math.Pow(fitParams[1], 2)) - y
  }
}
