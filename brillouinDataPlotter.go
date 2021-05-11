package main

import (
  "github.com/Arafatk/glot"
  "encoding/csv"
  "fmt"
  "os"
  "strconv"
)

func main() {


  label, freqFile, sigFile := readMeta()

  //data1 := getData(freqFile[1], sigFile[1])
  //data2 := getData(freqFile[2], sigFile[2])
  //data3 := getData(freqFile[3], sigFile[3])
  data4 := getData(freqFile[4], sigFile[4])
  //data5 := getData(freqFile[5], sigFile[5])
  //data6 := getData(freqFile[6], sigFile[6])
  //data7 := getData(freqFile[7], sigFile[7])
  //data8 := getData(freqFile[7], sigFile[7])
  //data9 := getData(freqFile[9], sigFile[9])


  plot(
    //data1, label[1],
    //data2, label[2],
    //data3, label[3],
    data4, label[4],
    //data5, label[5],
    //data6, label[6],
    //data7, label[7],
    //data8, label[7],
    //data9, label[9],
  )

}

func readMeta() ([]string, []string, []string) {

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

  var label, freq, sig []string

  for _, value := range meta {
    label = append(label, value[1])
    freq = append(freq, value[7])
    sig = append(sig, value[8])
  }

  return label, freq, sig
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
  set [][]float64, name string,
  //set2 [][]float64, name2 string,
  //set3 [][]float64, name3 string,
  //set4 [][]float64, name4 string,
  //set5 [][]float64, name5 string,
  //set6 [][]float64, name6 string,
  //set7 [][]float64, name7 string,
  //set8[][]float64, name8 string,
  //set9 [][]float64, name9 string,

  ) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("SBS in 1cm UHNA3 Fiber 50/50 Probe/LO -> 90/10")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (uV)")

  plot.AddPointGroup(name, "points", set)
  //plot.AddPointGroup(name2, "points", set2)
  //plot.AddPointGroup(name3, "points", set3)
  //plot.AddPointGroup(name4, "points", set4)
  //plot.AddPointGroup(name5, "points", set5)
  //plot.AddPointGroup(name6, "points", set6)
  //plot.AddPointGroup(name7, "points", set7)
  //plot.AddPointGroup(name8, "points", set8)
  //plot.AddPointGroup(name9, "points", set9)
}
