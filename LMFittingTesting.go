package main

import (
  "github.com/maorshutman/lm"
  "github.com/Arafatk/glot"
  "math"
  "fmt"
)

func main() {

  biggsNumJac := lm.NumJac{Func: biggsEXP6Func}

  biggsProb := lm.LMProblem{
	  Dim:        6,
 	  Size:       13,
 	  Func:       biggsEXP6Func,
 	  Jac:        biggsNumJac.Jac,
 	  InitParams: []float64{1, 2, 1, 1, 1, 1},
 	  Tau:        1e-6,
 	  Eps1:       1e-8,
 	  Eps2:       1e-8,
  }

  biggsResults, _ := lm.LM(biggsProb, &lm.Settings{Iterations: 100, ObjectiveTol: 1e-16})

  fmt.Println(biggsResults)

  fmt.Println(biggsResults.X)


  var z, y, yfit []float64

  for i := 0; i < 10; i++ {
    z = append(z, float64(i) / 10)
    y = append(y, math.Exp(-z[i]) - 5*math.Exp(-10*z[i]) + 3*math.Exp(-4*z[i]))

    yfit = append(yfit, biggsResults.X[2]*math.Exp(-biggsResults.X[0]*z[i]) - biggsResults.X[3]*math.Exp(-biggsResults.X[1]*z[i]) + biggsResults.X[5]*math.Exp(-biggsResults.X[4]*z[i]))
  }

  fmt.Println(z)
  fmt.Println(y)

  data := [][]float64{z, y}
  fit := [][]float64{z, yfit}

  plot, _ := glot.NewPlot(2, true, false)
  plot.AddPointGroup("data", "points", data)
  plot.AddPointGroup("fit", "lines", fit)

}

func biggsEXP6Func(dst, x []float64) {
	for i := 0; i < 13; i++ {
		z := float64(i) / 10
		y := math.Exp(-z) - 5*math.Exp(-10*z) + 3*math.Exp(-4*z)
		dst[i] = x[2]*math.Exp(-x[0]*z) - x[3]*math.Exp(-x[1]*z) + x[5]*math.Exp(-x[4]*z) - y
	}
}

/*
func powellFunc(dst, x []float64) {
	dst[0] = x[0]
	dst[1] = 10*x[0]/(x[0]+0.1) + 2*x[1]*x[1]
}
*/
