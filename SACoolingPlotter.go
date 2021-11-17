package main

import (
  "image/color"
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
  "gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
  "gonum.org/v1/plot/font"
	"gonum.org/v1/plot/vg/draw"
  //"gonum.org/v1/plot/vg/vgsvg"
  //"golang.org/x/image/font/gofont/goregular"
  //"github.com/golang/freetype/truetype"
)

func main() {

  label, file := readMeta()

  pras, pas, prs, ps := getAllData(file, label)
  prasLabel, pasLabel, prsLabel, psLabel := getAllLabels(label)

  toPlotRaw := []int{}
  if len(toPlotRaw) > 0 {
    plotRaw(
      toPlotRaw,
      pras, pas, prs, ps,
      prasLabel, pasLabel, prsLabel, psLabel,
    )
  }

  s, as := subtractBackground(pras, pas, prs, ps)

  toPlotSubtracted := []int{}
  if len(toPlotSubtracted) > 0 {
    plotSubtracted(toPlotSubtracted, s, as, prsLabel, prasLabel)
  }

  toPlotSubtractedTogether := []int{}
  if len(toPlotSubtractedTogether) > 0 {
    plotSubtractedTogether(toPlotSubtractedTogether, s, as, prsLabel, prasLabel)
  }

  toPlotSubtractedGrouped := []int{}
  if len(toPlotSubtractedGrouped) > 0 {
    plotSubtractedGrouped(toPlotSubtractedGrouped, s, as, prsLabel, prasLabel)
  }

  // gonum/plot
  numPlot := []int{0,1,2}
  if len(numPlot) > 0 {
    gonumPlot(numPlot, s, as, prsLabel, prasLabel)
  }

  // Lorentz fit better
  toFit := []int{}
  if len(toFit) > 0 {

    // Fit parameter guesses
    amp := 2.5
    wid := 0.2
    cen := 2.26

    // as
    fmt.Println("\nAnti-Stokes\n")
    var asFit [][]float64
    var asFits [][][]float64
    var asWidthLine [][]float64
    var asWidthLines [][][]float64
    var asfwhm []float64

    for _, set := range toFit {

      f := func(dst, guess []float64) {

        amp, wid, cen := guess[0], guess[1], guess[2]

        for i := 0; i < len(as[set][0]); i++ {
          x := as[set][0][i]
          y := as[set][1][i]
          dst[i] = amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)) - y
        }
      }

      jacobian := NumJac{Func: f}

      // Solve for fit
      toBeSolved := LMProblem{
    	  Dim:        3,
     	  Size:       len(as[set][0]),
     	  Func:       f,
     	  Jac:        jacobian.Jac,
     	  InitParams: []float64{amp, wid, cen},
     	  Tau:        1e-6,
     	  Eps1:       1e-8,
     	  Eps2:       1e-8,
      }

      results, _ := LM(toBeSolved, &Settings{Iterations: 100, ObjectiveTol: 1e-16})

      amp, wid, cen := results.X[0], math.Abs(results.X[1]), results.X[2]

      asfwhm = append(asfwhm, wid*2000)

      fmt.Printf("set %d | width: %.2f MHz | peak: %.6f nV | center: %.4f GHz\n", set, wid*2000, amp, cen)

      var asyFits []float64

      // Create function according to solved fit parameters
      for i := 0; i < len(as[set][0]); i++ {
          // (amp*wid^2/((x-cen)^2+wid^2))
          x := as[set][0][i]
          asyFits = append(asyFits, amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)))
      }

      // Width lines
      asWidthLine = [][]float64{{cen - wid, cen + wid},{amp/2, amp/2}}
      asWidthLines = append(asWidthLines, asWidthLine)

      asFit = [][]float64{as[set][0], asyFits}
      asFits = append(asFits, asFit)
    }

    // Plot fit
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Anti-Stokes")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (nV)")

    for _, set := range toFit {
      plot.AddPointGroup(strings.Trim(prasLabel[set], " pras"), "points", as[set])
      plot.AddPointGroup(strings.Trim(prasLabel[set], " pras") + " fit", "lines", asFits[set])
      plot.AddPointGroup(strconv.FormatFloat(asfwhm[set], 'f', 1, 64) + " MHz", "lines", asWidthLines[set])
    }

    // Width vs pump power
    asWidthPoints := [][]float64{{0, 80, 170},{asfwhm[0], asfwhm[1], asfwhm[2]}}

    // s
    fitStokes := []int{}
    if len(fitStokes) > 0 {

      fmt.Println("\nStokes\n")
      var sFit [][]float64
      var sFits [][][]float64
      var sWidthLine [][]float64
      var sWidthLines [][][]float64
      var sfwhm []float64

      for _, set := range toFit {

        f := func(dst, guess []float64) {

          amp, wid, cen := guess[0], guess[1], guess[2]

          for i := 0; i < len(s[set][0]); i++ {
            x := s[set][0][i]
            y := s[set][1][i]
            dst[i] = amp * math.Pow(wid, 2) / (math.Pow(x - cen, 2) + math.Pow(wid, 2)) - y
          }
        }

        jacobian := NumJac{Func: f}

        // Solve for fit
        toBeSolved := LMProblem{
      	  Dim:        3,
       	  Size:       len(s[set][0]),
       	  Func:       f,
       	  Jac:        jacobian.Jac,
       	  InitParams: []float64{amp, wid, cen},
       	  Tau:        1e-6,
       	  Eps1:       1e-8,
       	  Eps2:       1e-8,
        }

        results, _ := LM(toBeSolved, &Settings{Iterations: 100, ObjectiveTol: 1e-16})

        amp, wid, cen := results.X[0], math.Abs(results.X[1]), results.X[2]

        sfwhm = append(sfwhm, wid*2000)

        fmt.Printf("set %d | width: %.2f MHz | peak: %.6f nV | center: %.4f GHz\n", set, wid*2000, amp, cen)

        var syFits []float64

        // Create function according to solved fit parameters
        for i := 0; i < len(s[set][0]); i++ {
            // (amp*wid^2/((x-cen)^2+wid^2))
            x := s[set][0][i]
            syFits = append(syFits, results.X[0] * math.Pow(results.X[1], 2) / (math.Pow(x - results.X[2], 2) + math.Pow(results.X[1], 2)))
        }

        // Width lines
        sWidthLine = [][]float64{{cen - wid, cen + wid},{amp/2, amp/2}}
        sWidthLines = append(sWidthLines, sWidthLine)

        sFit = [][]float64{s[set][0], syFits}
        sFits = append(sFits, sFit)
      }
      fmt.Printf("\n")

      // Plot fit
      dimensions = 2
      persist = true
      debug = false
      plot, _ = glot.NewPlot(dimensions, persist, debug)

      plot.SetTitle("Stokes")
      plot.SetXLabel("Frequency (GHz)")
      plot.SetYLabel("Signal (nV)")

      for _, set := range toFit {
        //plot.AddPointGroup(strings.Trim(prsLabel[set], " prs"), "points", s[set])
        plot.AddPointGroup(strings.Trim(prsLabel[set], " prs"), "lines", sFits[set])
        plot.AddPointGroup(strconv.FormatFloat(sfwhm[set], 'f', 1, 64) + " MHz", "lines", sWidthLines[set])
      }

      // Width vs pump power
      //sWidthPoints := [][]float64{{0, 80, 170},{sfwhm[0], sfwhm[1], sfwhm[2]}}
    }

    // Plot width points
    dimensions = 2
    persist = true
    debug = false
    plot, _ = glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Pump Power vs Width")
    plot.SetXLabel("Pump Power (mW)")
    plot.SetYLabel("FWHM (MHz)")

    plot.AddPointGroup("Anti-Stokes", "circle", asWidthPoints)
    //plot.AddPointGroup("Stokes", "circle", sWidthPoints)

  }
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

func getAllData(fileNames, labels []string) ([][][]float64, [][][]float64, [][][]float64, [][][]float64) {

  var pras, pas, prs, ps [][][]float64

  // Assign data by checking csv name
  for _, fileName := range fileNames {
    if strings.Contains(fileName, "pras.") {
      pras = append(pras, getData(&fileName))
    } else if strings.Contains(fileName, "pas.") {
      pas = append(pas, getData(&fileName))
    } else if strings.Contains(fileName, "prs.") {
      prs = append(prs, getData(&fileName))
    } else if strings.Contains(fileName, "ps.") {
      ps = append(ps, getData(&fileName))
    }
  }

  return pras, pas, prs, ps
}

func getData(csvName *string) ([][]float64) {

  // Read
  f, err := os.Open(*csvName)
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

  for i := 1; i < 602; i++ {
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
    var nV []float64

    for _, dBm := range signal {
      nV = append(nV, 1000*math.Pow(10, 6)*math.Pow(10, dBm/10.))
    }

    return [][]float64{frequency, nV}
  } else if dataStr[1][3] == "  uV" {
    var nV []float64

    for _, uV := range signal {
      nV = append(nV, 1000*uV)
    }

    return [][]float64{frequency, nV}
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

  var prasLabel, pasLabel, prsLabel, psLabel []string

  // Assign labels by checking label
  for _, thisLabel := range label {
    if strings.Contains(thisLabel, "pras") {
      prasLabel = append(prasLabel, thisLabel)
    } else if strings.Contains(thisLabel, "pas") {
      pasLabel = append(pasLabel, thisLabel)
    } else if strings.Contains(thisLabel, "prs") {
      prsLabel = append(prsLabel, thisLabel)
    } else if strings.Contains(thisLabel, "ps") {
      psLabel = append(psLabel, thisLabel)
    }
  }

  return prasLabel, pasLabel, prsLabel, psLabel
}

func plotRaw(
  sets []int,
  pras, pas, prs, ps [][][]float64,
  prasLabel, pasLabel, prsLabel, psLabel []string,
  ) {

  for i := 0; i < len(sets); i++ {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Raw")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (uV)")

    plot.AddPointGroup(prasLabel[sets[i]], "points", pras[sets[i]])
    plot.AddPointGroup(pasLabel[sets[i]], "points", pas[sets[i]])
    plot.AddPointGroup(prsLabel[sets[i]], "points", prs[sets[i]])
    plot.AddPointGroup(psLabel[sets[i]], "points", ps[sets[i]])
  }
}

func subtractBackground(pras, pas, prs, ps [][][]float64) ([][][]float64, [][][]float64) {

  var s, as [][][]float64

  for i := 0; i < len(pras); i++ {
    s = append(s, subtract(ps[i], prs[i]))
    as = append(as, subtract(pas[i], pras[i]))
  }

  return s, as
}

func subtract(b, s [][]float64) ([][]float64) {

  var sum float64
  n := 100

  for i := 0; i < n; i++ {
    sum += b[1][i] - s[1][i]
  }

  for i := 0; i < len(b[0]); i++ {
    s[1][i] = s[1][i] - b[1][i] + sum/float64(n)
  }

  return s
}

func plotSubtracted(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  for _, set := range sets {
    dimensions := 2
    persist := true
    debug := false
    plot, _ := glot.NewPlot(dimensions, persist, debug)

    plot.SetTitle("Background Subtracted")
    plot.SetXLabel("Frequency (GHz)")
    plot.SetYLabel("Signal (nV)")

    plot.AddPointGroup(strings.Trim(sLabel[set], " prs") + " s", "points", s[set])
    plot.AddPointGroup(strings.Trim(asLabel[set], " pras") + " as", "points", as[set])
  }
}

func plotSubtractedTogether(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (nV)")

  for _, set := range sets {
    plot.AddPointGroup(strings.Trim(sLabel[set], " prs") + " s", "points", s[set])
    plot.AddPointGroup(strings.Trim(asLabel[set], " pras") + " as", "points", as[set])
  }
}

func plotSubtractedGrouped(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  // s
  dimensions := 2
  persist := true
  debug := false
  plot, _ := glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (nV)")

  for _, set := range sets {
    plot.AddPointGroup(strings.Trim(sLabel[set], " prs") + " s", "points", s[set])
  }

  // as
  dimensions = 2
  persist = true
  debug = false
  plot, _ = glot.NewPlot(dimensions, persist, debug)

  plot.SetTitle("Background Subtracted")
  plot.SetXLabel("Frequency (GHz)")
  plot.SetYLabel("Signal (nV)")

  for _, set := range sets {
    plot.AddPointGroup(strings.Trim(asLabel[set], " pras") + " as", "points", as[set])
  }
}

func gonumPlot(sets []int, s, as [][][]float64, sLabel, asLabel []string) {

  /*goFont, err := truetype.Parse(goregular.TTF)
  if err != nil {
    panic(err)
  }*/

  // as
  p := plot.New()
  p.Title.Text = "Anti-Stokes"
  p.Title.TextStyle.Font.Typeface = "liberation"
  p.Title.TextStyle.Font.Variant = "Sans"
  p.Title.TextStyle.Font.Size = 34
  p.Title.Padding = font.Length(50)

  p.X.Label.Text = "Frequency (GHz)"
  p.X.Label.TextStyle.Font.Variant = "Sans"
  p.X.Label.TextStyle.Font.Size = 24
  p.X.Label.Padding = font.Length(20)
  p.X.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.LineStyle.Width = vg.Points(1.5)
  p.X.Tick.Label.Font.Size = 24
  p.X.Tick.Label.Font.Variant = "Sans"
  p.X.Tick.Marker = plot.ConstantTicks([]plot.Tick{
  		{Value: 2.0, Label: "2"},
      {Value: 2.05, Label: ""},
      {Value: 2.1, Label: "2.1"},
      {Value: 2.15, Label: ""},
      {Value: 2.2, Label: "2.2"},
      {Value: 2.25, Label: ""},
      {Value: 2.3, Label: "2.3"},
      {Value: 2.35, Label: ""},
      {Value: 2.4, Label: "2.4"},
      {Value: 2.45, Label: ""},
      {Value: 2.5, Label: "2.5"},
  	})
  p.X.Padding = vg.Points(25)

  p.Y.Label.Text = "Signal (nV)"
  p.Y.Label.TextStyle.Font.Variant = "Sans"
  p.Y.Label.TextStyle.Font.Size = 24
  p.Y.Label.Padding = font.Length(20)
  p.Y.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.LineStyle.Width = vg.Points(1.5)
  p.Y.Tick.Label.Font.Size = 24
  p.Y.Tick.Label.Font.Variant = "Sans"
  p.Y.Tick.Marker = plot.ConstantTicks([]plot.Tick{
  		{Value: -1, Label: "-1"},
      {Value: -.75, Label: ""},
      {Value: -.5, Label: "-.5"},
      {Value: -.25, Label: ""},
      {Value: 0, Label: "0"},
      {Value: .25, Label: ""},
      {Value: .5, Label: ".5"},
      {Value: .75, Label: ""},
      {Value: 1, Label: "1"},
      {Value: 1.25, Label: ""},
      {Value: 1.5, Label: "1.5"},
      {Value: 1.75, Label: ""},
      {Value: 2, Label: "2"},
      {Value: 2.25, Label: ""},
      {Value: 2.5, Label: "2.5"},
      {Value: 2.75, Label: ""},
      {Value: 3, Label: "3"},

  	})
  p.Y.Padding = vg.Points(25)

  p.Legend.TextStyle.Font.Size = 24
  p.Legend.TextStyle.Font.Variant = "Sans"
  p.Legend.Top = true
  p.Legend.XOffs = vg.Points(-50)
  p.Legend.YOffs = vg.Points(-50)
  p.Legend.Padding = vg.Points(10)
  p.Legend.ThumbnailWidth = vg.Points(50)

  setColors := make([]color.RGBA, len(sets))
  setColors[0] = color.RGBA{R: 31, G: 249, B: 155, A: 255}
  setColors[1] = color.RGBA{R: 255, G: 122, B: 180, A: 255}
  setColors[2] = color.RGBA{R: 122, G: 156, B: 255, A: 255}


  for _, set := range sets {

    asPts := buildData(as[set])

  	// Make a scatter plotter and set its style.
  	plotSet, err := plotter.NewScatter(asPts)
  	if err != nil {
  		panic(err)
  	}

  	plotSet.GlyphStyle.Color = setColors[set]
    plotSet.GlyphStyle.Radius = vg.Points(3)
    plotSet.Shape = draw.CircleGlyph{}

    p.Add(plotSet)
    p.Legend.Add(strings.Trim(asLabel[set], " pras"), plotSet)
  }

  savePlotAs := "Anti-Stokes Background Subtracted"
  // Save the plot to a PNG file.
  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".png"); err != nil {
  	panic(err)
  }

  if err := p.Save(15*vg.Inch, 15*vg.Inch, savePlotAs+".svg"); err != nil {
    panic(err)
  }
}

func buildData(data [][]float64) (plotter.XYs) {

  xy := make(plotter.XYs, len(data[0]))

  for i := range xy {
    xy[i].X = data[0][i]
    xy[i].Y = data[1][i]
  }

  return xy
}

func normalizeFit(fit []float64) ([]float64) {

  var shift float64 = (fit[0] + fit[599])/2

  for i := 0; i < 600; i++ {
    fit[i] = fit[i] - shift
  }
  return fit
}

//----------------------------------------------------------------------------\\

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
	Func func(dst, guess []float64)
}

func (nj *NumJac) Jac(dst *mat.Dense, guess []float64) {
	fd.Jacobian(dst, nj.Func, guess, &fd.JacobianSettings{
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
