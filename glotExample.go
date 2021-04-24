package main

import "github.com/Arafatk/glot"

func main() {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    circles := [][]float64{{7, 3, 13, 5.6, 11.1}, {12, 13, 11, 1,  7}}
        plot.AddPointGroup("Circles", "circle", circles)
    points := [][]float64{{3, 6, 1, 5.5, 12}, {7.6, 8, 9, 3, 9}}
        plot.AddPointGroup("Points", "points", points)
    lines := [][]float64{{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {2, 9, 2, 4.5, 7, 4, 3, 6.7, 8, 2.3}}
        plot.AddPointGroup("Lines", "lines", lines)
    impulses := [][]float64{{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {9, 2, 4.5, 7, 4, 3, 6.7, 8, 2.3, 5}}
        plot.AddPointGroup("Impulses", "impulses", impulses)
    points2 := [][]float64{{4, 5, 6, 7, 8, 9}, {9, 2, 4.5, 7, 4, 3}}
        plot.AddPointGroup("Points2", "points", points2)
    points3 := [][]float64{{5, 6, 7, 8, 9}, {9, 2, 4.5, 7, 4}}
        plot.AddPointGroup("Points3", "points", points3)

    plot.SetTitle("Example Plot")
    plot.SetXLabel("X-Axis")
    plot.SetYLabel("Y-Axis")
    plot.SetXrange(-2, 18)
    plot.SetYrange(-2, 18)
    plot.SavePlot("Scatterplot.png")
}
