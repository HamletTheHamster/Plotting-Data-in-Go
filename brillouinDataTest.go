package main

import (
  "github.com/Arafatk/glot"
  "encoding/csv"
  "fmt"
  "os"
  "strconv"
)

func main() {
  data100 := getData(
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/100%/100% frequency.csv",
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/100%/100% signal.csv",
  )
  data30 := getData(
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/30%/30% frequency.csv",
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/30%/30% signal.csv",
  )
  data15 := getData(
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/15%/15% frequency.csv",
    "../../../Data/2 - 1cm UHNA3 v2/2021-3-5/15%/15% signal.csv",
  )
  plot(
    data100, "100%",
    data30, "30%",
    data15, "15%",
  )
}

func getData(frequencyCSV string, signalCSV string) ([][]float64) {

  // Read
  file, err := os.Open(frequencyCSV)
  if err != nil {
      fmt.Println(err)
  }

  reader := csv.NewReader(file)
  frequencyStr, _ := reader.ReadAll()

  file, err = os.Open(signalCSV)
  if err != nil {
      fmt.Println(err)
  }

  reader = csv.NewReader(file)
  signalStr, _ := reader.ReadAll()

  // Transpose
  var frequencyStrT, signalStrT []string

  for i := 0; i < 300; i++ {
    frequencyStrT = append(frequencyStrT, frequencyStr[i][0])
    signalStrT = append(signalStrT, signalStr[i][0])
  }

  // Convert
  var frequency, signal []float64

  for _, freqElem := range frequencyStrT {
    freqValue, err := strconv.ParseFloat(freqElem, 64)
    if err == nil {
      frequency = append(frequency, freqValue/1e9)
    }
  }

  for _, sigElem := range signalStrT {
    sigValue, err := strconv.ParseFloat(sigElem, 64)
    if err == nil {
      signal = append(signal, sigValue*1e6)
    }
  }

  return [][]float64{frequency,signal}
}

func plot(
  set1 [][]float64, name1 string,
  set2 [][]float64, name2 string,
  set3 [][]float64, name3 string,
  ) {
  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("SBS in 1cm UHNA3 Fiber")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uW)")
  plot.AddPointGroup(name1+"%", "points", set1)
  plot.AddPointGroup(name2+"%", "points", set2)
  plot.AddPointGroup(name3+"%", "points", set3)
}
