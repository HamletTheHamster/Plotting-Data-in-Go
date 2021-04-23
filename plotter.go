package main

import (
    "math/rand"

    "gonum.org/v1/plot"
    "gonum.org/v1/plot/plotter"
    "gonum.org/v1/plot/vg"
)

func main() {
    //make data
    var values plotter.Values
    for i := 0; i < 1000; i++ {
        values = append(values, rand.NormFloat64())
    }
    //boxPlot(values)
    //barPlot(values[:4])
    histPlot(values)

}

func histPlot(values plotter.Values) {
    p, err := plot.New()
    if err != nil {
        panic(err)
    }
    p.Title.Text = "histogram plot"

    hist, err := plotter.NewHist(values, 20)
    if err != nil {
        panic(err)
    }
    p.Add(hist)

    if err := p.Save(3*vg.Inch, 3*vg.Inch, "hist.png"); err != nil {
        panic(err)
    }
}
