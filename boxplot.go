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
    boxPlot(values[:4])

}

func boxPlot(values plotter.Values) {
    p := plot.New()
    p.Title.Text = "box plot"

    box, err := plotter.NewBoxPlot(vg.Length(15), 0.0, values)
    if err != nil {
        panic(err)
    }
    p.Add(box)

    if err := p.Save(9*vg.Inch, 9*vg.Inch, "box.png"); err != nil {
        panic(err)
    }
}
