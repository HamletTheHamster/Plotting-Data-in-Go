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
    barPlot(values[:4])

}

func barPlot(values plotter.Values) {
    p := plot.New()
    p.Title.Text = "bar plot"

    bar, err := plotter.NewBarChart(values, 15)
    if err != nil {
        panic(err)
    }
    p.Add(bar)

    if err := p.Save(9*vg.Inch, 9*vg.Inch, "bar.png"); err != nil {
        panic(err)
    }
}
